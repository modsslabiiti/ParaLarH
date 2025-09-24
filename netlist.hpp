#ifndef NETLIST_H
#define NETLIST_H

#include "placement.hpp"
#include "pugixml.hpp"
#include "basic_types.hpp"
#include "router.hpp"
#include "arch.hpp"
#include <boost/algorithm/string/trim.hpp>		/* "/usr/include/boost......' */

using namespace boost;
using namespace pugi;

struct Net {
	int index;
	string name;
	vector<Point> points;
	int *clusterid;
	double **cdata;
	int num_clusters;
	list<GridGraph::edge_descriptor> route_edges;
};

/*vector<Net> read_netlist(const string &filename);*/
void read_netlist_xml(const string &filename);

namespace boost {
	enum vertex_tnode_t { vertex_tnode };

	BOOST_INSTALL_PROPERTY(vertex, tnode);
}

struct Pin {
	const struct Port *port;
	int index;
	string driver;
	int tnode;
};

struct Port {
	const struct Block *block;
	string name;
	vector<std::unique_ptr<Pin>> pins;
	ptrdiff_t file_offset;
};

struct Block {
	string name;
	string instance;
	int index;
	Point position;
	Pb *pb;
	Mode *mode;
	const Block *parent;
	map<string, std::unique_ptr<Port>> input_ports;
	map<string, std::unique_ptr<Port>> output_ports;
	map<string, std::unique_ptr<Port>> clock_ports;
	map<string, std::unique_ptr<Block>> children;
};

struct NewNet {
	string name;
	const Pin *source;
	vector<const Pin *> sinks;
};

struct timing_node {
	const Pin *pin;
	struct timing_node *prev;
	float arrival_time;
	float required_time;
};

class Netlist {
	public:
		map<string, std::unique_ptr<NewNet>> nets;

		void read_port(const xml_node &block_node, Block *b)
		{
			vector<string> types { "inputs", "outputs", "clocks" };
			for (auto t : types) {
				auto node = block_node.child(t.c_str());
				for (auto port_node = node.child("port"); port_node; port_node = port_node.next_sibling("port")) {
					std::unique_ptr<Port> port { new Port() };
					string tmp = port_node.text().get();
					boost::trim_if(tmp, boost::is_any_of("\t "));
					vector<string> tokens;
					boost::split(tokens, tmp, boost::is_any_of(" "), boost::token_compress_on);

					int index = 0;
					for (auto token : tokens) {
						std::unique_ptr<Pin> pin { new Pin() };
						pin->driver = token;
						if (token != "open") {
							pin->tnode = add_vertex(tgraph);
							tgraph[pin->tnode].pin = pin.get();
						} else {
							pin->tnode = graph_traits<__decltype(tgraph)>::null_vertex();
						}
						pin->index = index++;
						pin->port = port.get();
/*						put(vertex_tnode, tgraph, pin->tnode, pin.get());*/
/*						printf("%s ", pin->driver.c_str());*/
						port->pins.emplace_back(std::move(pin));
					}
/*					printf("\n");*/

					port->name = port_node.attribute("name").as_string();
					port->block = b;
					port->file_offset = port_node.offset_debug();

					if (t == "inputs") {
						assert(b->input_ports.find(port->name) == b->input_ports.end());
/*						printf("adding input port %s\n", port->name.c_str());*/
						b->input_ports[port->name] = std::move(port);
					} else if (t == "outputs") {
						assert(b->output_ports.find(port->name) == b->output_ports.end());
/*						printf("adding output port %s\n", port->name.c_str());*/
						b->output_ports[port->name] = std::move(port);
					} else if (t == "clocks") {
						assert(b->clock_ports.find(port->name) == b->clock_ports.end());
/*						printf("adding output port %s\n", port->name.c_str());*/
						b->clock_ports[port->name] = std::move(port);
					} else {
						assert(false);
					}
				}
			}
		}

		void resolve_input_port(const Block *b)
		{
			static boost::regex pattern(R"d(((\w+)(?:\[(\d+)\]){0,1})\.((\w+)(?:\[(\d+)\]))->(.*))d");

			for (const auto &port : b->input_ports) {
				for (const auto &pin : port.second->pins) {
					if (pin->driver != "open") {
						boost::smatch match;
						if (boost::regex_search(pin->driver, match, pattern)) {
//							for (auto &m : match) {
//								printf("%s \n", m.str().c_str());
//							}
							if (match.size() == 8) { 
								const Block *driver_block = nullptr;
								map<string, std::unique_ptr<Port>>::const_iterator driver_port;
								if (match.str(3).empty()) { /* from parent */
/*									printf("from parent\n");*/
									driver_block = b->parent;
									assert(driver_block->instance.find(match.str(1)) == 0 && match.str(1) == match.str(2));
									driver_port = driver_block->input_ports.find(match.str(5));
									assert(driver_port != driver_block->input_ports.end());
								} else { /*from sibling*/
/*									printf("from sibling\n");*/
									string child_instance = match.str(1);
									boost::trim(child_instance);
									auto child = b->parent->children.find(child_instance);
									assert(child != b->parent->children.end());
									driver_block = child->second.get();
									assert(driver_block->instance == match.str(1));
									driver_port = driver_block->output_ports.find(match.str(5));
									assert(driver_port != driver_block->output_ports.end());
								}
/*								printf("looking for %s from %s %s\n", match.str(5).c_str(), driver_block->name.c_str(), driver_block->instance.c_str());*/

								auto driver_pin = driver_port->second->pins[boost::lexical_cast<int>(match.str(6).c_str())].get();
								auto test = edge(driver_pin->tnode, pin->tnode, tgraph);
								assert(!test.second);
								auto res = add_edge(driver_pin->tnode, pin->tnode, tgraph);
								assert(res.second == true);

								float d = 0;
								int num_matches = 0;
								bool found = false;
								for (const auto &interconnect : b->parent->mode->interconnects) {
									if (interconnect->name == match.str(7)) {
										int num_inner_matches = 0;
										for (const auto &delay : interconnect->delays) {
											int in_instance_index;
											if (match.str(3) == "") {
												in_instance_index = 0;
											} else {
												in_instance_index = boost::lexical_cast<int>(match.str(3));
											}
											if (delay->matches(match.str(2), in_instance_index, match.str(5), boost::lexical_cast<int>(match.str(6)), b->instance.substr(0, b->instance.find('[')), b->index, port.second->name, pin->index)) {
												assert(delay->type == DELAY_CONSTANT);

												d = delay->matrix[0][0];
												//printf("%s[%d].%s[%d] -> %s[%d].%s[%d] Delay: %g\n", match.str(2).c_str(), in_instance_index, match.str(5).c_str(), boost::lexical_cast<int>(match.str(6)), b->instance.substr(0, b->instance.find('[')).c_str(), b->index, port.second->name.c_str(), pin->index, d);
/*												printf("%s.%s -> %s.%s Delay: %g\n", match.str(2).c_str(), match.str(5).c_str(), b->instance.substr(0, b->instance.find('[')).c_str(), port.second->name.c_str(), d);*/

												if (found) {
													assert(false);
												} else {
													found = true;
												}
												num_inner_matches++;
											}
										}
										assert(num_inner_matches <= 1);
										num_matches++;
									}
								}
								assert(num_matches > 0);

								put(edge_weight, tgraph, res.first, d);
							} else {
								assert(false);
							}
						} else {
							/*top level input ports*/
							assert(b->parent && !b->parent->parent);
							auto net = nets.find(pin->driver);
							if (net == nets.end()) {
								std::unique_ptr<NewNet> new_net { new NewNet() };
								new_net->name = pin->driver;
								new_net->sinks.emplace_back(pin.get());

								nets[pin->driver] = std::move(new_net);
							} else {
								net->second->sinks.emplace_back(pin.get());
							}
						}
					}
				}
			}
		}

		void resolve_top_level_output_port(const Block *b)
		{
			assert(b->parent && !b->parent->parent);

			for (const auto &port : b->output_ports) {
				for (const auto &pin : port.second->pins) {
					if (pin->driver != "open") {
						const Pin *cur_pin = pin.get();
						while (!cur_pin->port->block->children.empty()) {
							assert(in_degree(cur_pin->tnode, tgraph) == 1);
							auto edges = in_edges(cur_pin->tnode, tgraph);
							auto e = *edges.first;
							assert(target(e, tgraph) == cur_pin->tnode);
/*							cur_pin = get(vertex_tnode, tgraph, source(e, tgraph));*/
							cur_pin = tgraph[source(e, tgraph)].pin;
						}
						assert(cur_pin->port->block->instance.find("lut[") == 0 || cur_pin->port->block->instance.find("ff[") == 0 || cur_pin->port->block->instance.find("inpad[") == 0);
						auto net = nets.find(cur_pin->driver);
						if (net == nets.end() && cur_pin->driver != "pclk" && cur_pin->driver != "clock") {
							printf("can't find net %s\n", cur_pin->driver.c_str());
						}
						assert(net != nets.end() || cur_pin->driver == "pclk" || cur_pin->driver == "clock");
						net->second->source = pin.get();
					}
				}
			}
		}

		void resolve_output_port(const Block *b)
		{
			static boost::regex pattern(R"d(((\w+)(?:\[(\d+)\]){0,1})\.((\w+)(?:\[(\d+)\]))->(.*))d");

			for (const auto &port : b->output_ports) {
				for (const auto &pin : port.second->pins) {
					if (pin->driver != "open") {
						boost::smatch match;
						if (boost::regex_search(pin->driver, match, pattern)) {
/*							for (auto &m : match) {*/
/*								printf("%s \n", m.str().c_str());*/
/*							}	*/
							if (match.size() == 8) { 
								map<string, std::unique_ptr<Port>>::const_iterator driver_port;
								if (match.str(1) == b->instance) { /* from current */
/*									printf("from parent\n");*/
									auto driver_block = b;
									assert(driver_block->instance.find(match.str(1)) == 0);// && match.str(1) == match.str(2));
									driver_port = driver_block->input_ports.find(match.str(5));
									assert(driver_port != driver_block->input_ports.end());
								} else { /*from children */
/*									printf("from child\n");*/
									string child_instance = match.str(1);
									boost::trim(child_instance);
									auto child = b->children.find(child_instance);
									if (child == b->children.end()) {
										printf("offset: %lu\n", pin->port->file_offset);
										printf("children of %s %s\n", b->name.c_str(), b->instance.c_str());
										for (const auto &c : b->children) {
											printf("%s\n", c.first.c_str());
										}
										printf("can't find %s\n", child_instance.c_str());
									}
									assert(child != b->children.end());
									auto driver_block = child->second.get();
									if (driver_block->instance != match.str(1)) {
										printf("Name: %s Instance: %s\n", pin->port->block->name.c_str(), pin->port->block->instance.c_str());
										printf("Driver: %s Found block instance: %s Expected block instance: %s\n", pin->driver.c_str(), driver_block->instance.c_str(), match.str(1).c_str());
									}
									assert(driver_block->instance == match.str(1));
									driver_port = driver_block->output_ports.find(match.str(5));
									assert(driver_port != driver_block->output_ports.end());
								}
/*								printf("looking for %s from %s %s\n", match.str(5).c_str(), driver_block->name.c_str(), driver_block->instance.c_str());*/

								auto driver_pin = driver_port->second->pins[boost::lexical_cast<int>(match.str(6).c_str())].get();
								auto test = edge(driver_pin->tnode, pin->tnode, tgraph);
								assert(!test.second);
								auto res = add_edge(driver_pin->tnode, pin->tnode, tgraph);
								assert(res.second == true);

								float d = 0;
								int num_matches = 0;
								bool found = false;
								for (const auto &interconnect : b->mode->interconnects) {
									if (interconnect->name == match.str(7)) {
										int num_inner_matches = 0;
										for (const auto &delay : interconnect->delays) {
											int in_instance_index;
											if (match.str(1) == b->instance) {
												in_instance_index = 0;
											} else {
												in_instance_index = boost::lexical_cast<int>(match.str(3));
											}
											if (delay->matches(match.str(2), in_instance_index, match.str(5), boost::lexical_cast<int>(match.str(6)), b->instance.substr(0, b->instance.find('[')), b->index, port.second->name, pin->index)
											   ) {
												assert(delay->type == DELAY_CONSTANT);

												d = delay->matrix[0][0];
												//printf("%s[%d].%s[%d] -> %s[%d].%s[%d] Delay: %g\n", match.str(2).c_str(), in_instance_index, match.str(5).c_str(), boost::lexical_cast<int>(match.str(6)), b->instance.substr(0, b->instance.find('[')).c_str(), b->index, port.second->name.c_str(), pin->index, d);
/*												printf("%s.%s -> %s.%s Delay: %g\n", match.str(2).c_str(), match.str(5).c_str(), b->instance.substr(0, b->instance.find('[')).c_str(), port.second->name.c_str(), d);*/
												if (found) {
													assert(false);
												} else {
													found = true;
												}
												num_inner_matches++;
											}
										}
										assert(num_inner_matches <= 1);
										num_matches++;
									}
								}
								
								assert(num_matches > 0);
								put(edge_weight, tgraph, res.first, d);
							} else {
								assert(false);
							}
						} else {
							assert(b->children.empty());
/*							printf("Net driver: %s\n", pin->driver.c_str());*/
						}
					}
				}
			}
		}

		void check(const Block *b)
		{
			if (b->instance.find("ff[") == 0) {
				assert(b->children.empty());
/*				assert(false);*/
/*				for (*/
				auto iter = b->output_ports.find("Q");
				if (iter != b->output_ports.end()) {
					auto output_pin = iter->second->pins[0].get();
					for (const auto &port : b->input_ports) {
						for (const auto &pin : port.second->pins) {
							if (pin->driver != "open") {
								auto test = edge(pin->tnode, output_pin->tnode, tgraph);
								assert(!test.second);
								test = edge(output_pin->tnode, pin->tnode, tgraph);
								assert(!test.second);
							}
						}
					}
				} else {
					assert(b->name == "open");
				}
			}
			for (const auto &child : b->children) {
				check(child.second.get());
			}
		}

		void connect_primitives(const Block *b)
		{
			if (b->instance.find("lut[") == 0) {
				assert(b->children.empty());
				assert(b->output_ports.size() == 1 && b->output_ports.find("out")->second->pins.size() == 1);
				auto output_pin = b->output_ports.find("out")->second->pins[0].get();
				
				for (const auto &port : b->input_ports) {
					for (const auto &pin : port.second->pins) {
						if (pin->driver != "open") {
							auto test = edge(pin->tnode, output_pin->tnode, tgraph);
							assert(!test.second);
							auto res = add_edge(pin->tnode, output_pin->tnode, tgraph);
							assert(res.second == true);
							put(edge_weight, tgraph, res.first, 235e-12);
						}
					}
				}
			}
		}

		void resolve(const Block *b)
		{
			/* don't resolve for root block */
			if (b->parent) {
				resolve_input_port(b);
				resolve_output_port(b);
				connect_primitives(b);
			}
			for (const auto &child : b->children) {
				resolve(child.second.get());
			}
		}

		std::unique_ptr<Block> read_block(const xml_node &block_node, const Block *parent, const Placement &placement, const map<string, struct Pb *> &pbs)
		{
			std::unique_ptr<Block> b { new Block() };
			xml_attribute instance = block_node.attribute("instance");
			xml_attribute name = block_node.attribute("name");
			auto mode_name = block_node.attribute("mode");
			b->instance = instance.as_string();
			int start = b->instance.find('[');
			int end = b->instance.find(']');
			b->index = boost::lexical_cast<int>(b->instance.substr(start+1, end-start+1-2));
			b->name = name.as_string();
			b->parent = parent;
			if (b->parent && !b->parent->parent) {
				auto pos = placement.getBlockPosition(b->name);
				b->position.x = pos.first;
				b->position.y = pos.second;
				assert(b->position.x >= 0 && b->position.y >= 0);
			}

			read_port(block_node, b.get());

/*			printf("name: %s instance: %s\n", name.as_string(), instance.as_string());*/

			if (mode_name) {
				auto tmp = b->instance.substr(0, b->instance.find('['));
				auto iter = pbs.find(tmp);
				assert(iter != pbs.end());
				b->pb = iter->second;
				auto mode_iter = b->pb->modes.find(mode_name.as_string());
				assert(mode_iter != b->pb->modes.end());
				b->mode = mode_iter->second;
				assert(b->mode->name == mode_name.as_string());
			} else {
				b->pb = nullptr;
				b->mode = nullptr;
				assert(parent == NULL || block_node.child("block").empty());
			}

			for (auto child_node = block_node.child("block"); child_node; child_node = child_node.next_sibling("block")) {
				std::unique_ptr<Block> child_block;
				if (b->mode) {
					child_block = read_block(child_node, b.get(), placement, b->mode->pbs);
				} else {
					child_block = read_block(child_node, b.get(), placement, pbs);
				}

				assert(b->children.find(child_block->instance) == b->children.end());
				string instance = child_block->instance;
				b->children[instance] = std::move(child_block);
			}
			return b;
		}

		void load(const string &filename, const map<string, struct Pb *> &pbs)
		{
			xml_document doc;
			xml_parse_result result = doc.load_file((filename+".net").c_str());
			assert(result.status == status_ok);
			xml_node root_node = doc.document_element();

			Placement placement(filename + ".place");

			grid_size = placement.getGridSize();
			root_block = read_block(root_node, nullptr, placement, pbs);

			resolve(root_block.get());

			for (auto &child : root_block->children) {
				resolve_top_level_output_port(child.second.get());
			}

			for (const auto &net : nets) {
				for (const auto &sink : net.second->sinks) {
					auto test = edge(net.second->source->tnode, sink->tnode, tgraph);
					assert(test.second == false);
					add_edge(net.second->source->tnode, sink->tnode	, tgraph);
				}
			}

			for (auto &child : root_block->children) {
				check(child.second.get());
			}

			assert(root_block->name == filename+".net" && root_block->instance == "FPGA_packed_netlist[0]");
		}


/*		typedef property<vertex_tnode_t, const Pin *,*/
/*		       	property<vertex_bundle_t, timing_node>> TimingGraphVertexProperty;*/
		
/*		typedef property<vertex_bundle_t, int> TimingGraphVertexProperty;*/

		boost::adjacency_list<vecS, vecS, bidirectionalS, timing_node, property<edge_weight_t, float>> tgraph;
		std::unique_ptr<Block> root_block;
		pair<int, int> grid_size;

};

#endif
