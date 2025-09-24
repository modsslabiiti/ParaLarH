#include "pch.hpp"
#include "pugixml.hpp"
#include <stdarg.h>
#include "arch.hpp"

void print_tabs(int num_tabs, const char *str, ...)
{
/*	for (int i = 0; i < num_tabs; ++i) {*/
/*		printf("\t");*/
/*	}*/
/*	va_list args;*/
/*	va_start(args, str);*/
/*	vprintf(str, args);*/
/*	va_end(args);*/
}

bool Delay::matches(const string &in_instance, int in_instance_index, const string &in_port, int in_pin, const string &out_instance, int out_instance_index, const string &out_port, int out_pin) const
{
/*	printf("%s %d %s %d, %s %d %s %d\n", in_instance.c_str(), in_instance_index, in_port.c_str(), in_pin, out_instance.c_str(), out_instance_index, out_port.c_str(), out_pin);*/
/*	printf("%s %d:%d %s %d:%d, %s %d:%d %s %d:%d\n",*/
/*			this->in_port.instance.c_str(), this->in_port.instance_low, this->in_port.instance_high,*/
/*			this->in_port.port.c_str(), this->in_port.port_low, this->in_port.port_high,*/
/*			this->out_port.instance.c_str(), this->out_port.instance_low, this->out_port.instance_high,*/
/*						this->out_port.port.c_str(), this->out_port.port_low, this->out_port.port_high);*/
	return in_instance == this->in_port.instance && in_instance_index >= this->in_port.instance_low && in_instance_index <= this->in_port.instance_high &&
		in_port == this->in_port.port && in_pin >= this->in_port.port_low && in_pin <= this->in_port.port_high &&
		out_instance == this->out_port.instance && out_instance_index >= this->out_port.instance_low && out_instance_index <= this->out_port.instance_high &&
		out_port == this->out_port.port && out_pin >= this->out_port.port_low && out_pin <= this->out_port.port_high;
}

void load_low_high(const string &low_str, const string &high_str, int *low, int *high, bool instance)
{
	if (low_str != "") {
		*low = boost::lexical_cast<int>(low_str);
	} else {
		*low = 0;
	}

	if (high_str != "") {
		*high = boost::lexical_cast<int>(high_str);
	} else {
		if (instance) {
			*high = 0;
		} else {
			*high = INT_MAX;
		}
	}
	
	if (*high < *low) {
		int tmp = *high;
		*high = *low;
		*low = tmp;
	}
}

void parse_delay_input_output(const string &in_port, const string &out_port, Delay *delay)
{
	static boost::regex pattern(R"d(((\w+)(\[(\d+):(\d+)\]){0,1})\.(\w+)(\[(\d+):(\d+)\]){0,1})d");
	bool matches = false;

	boost::smatch in_match;
	if (boost::regex_search(in_port, in_match, pattern)) {
/*		printf("matches: %d\n", in_match.size());*/
		delay->in_port.instance = in_match.str(2);
		load_low_high(in_match.str(4), in_match.str(5), &delay->in_port.instance_low, &delay->in_port.instance_high, false);
		delay->in_port.port = in_match.str(6);
		load_low_high(in_match.str(8), in_match.str(9), &delay->in_port.port_low, &delay->in_port.port_high, false);
/*		printf("%s[%d:%d].%s[%d:%d]\n", delay->in_port.instance.c_str(), delay->in_port.instance_low, delay->in_port.instance_high,delay->in_port.port.c_str(), delay->in_port.port_low, delay->in_port.port_high);*/
/*		for (auto m : in_match) {*/
/*			printf("%s\n", m.str().c_str());*/
/*		}*/
	} else {
		assert(false);
	}

	boost::smatch out_match;
	if (boost::regex_search(out_port, out_match, pattern)) {
	/*		printf("matches: %d\n", in_match.size());*/
			delay->out_port.instance = out_match.str(2);
			load_low_high(out_match.str(4), out_match.str(5), &delay->out_port.instance_low, &delay->out_port.instance_high, false);
			delay->out_port.port = out_match.str(6);
			load_low_high(out_match.str(8), out_match.str(9), &delay->out_port.port_low, &delay->out_port.port_high, false);
/*			printf("%s[%d:%d].%s[%d:%d]\n", delay->out_port.instance.c_str(), delay->out_port.instance_low, delay->out_port.instance_high,delay->out_port.port.c_str(), delay->out_port.port_low, delay->out_port.port_high);*/
	/*		for (auto m : out_match) {*/
	/*			printf("%s\n", m.str().c_str());*/
	/*		}*/
		} else {
			assert(false);
		}
/*	printf("\n");*/

/*	boost::smatch output_match;*/
/*	if (matches && boost::regex_search(output, output_match, pattern)) {*/
/**/
/*	} else {*/
/*		matches = false;*/
/*	}*/
/*	return matches;*/
}

void Arch::parse_delays(const xml_node &interconnect_node, Interconnect *interconnect, int level)
{
	for (auto constant_node = interconnect_node.child("delay_constant"); constant_node; constant_node = constant_node.next_sibling("delay_constant")) {
		Delay *delay = new Delay();
		delay->type = DELAY_CONSTANT;
		delay->matrix = (float **)malloc(sizeof(float *));
		delay->matrix[0] = (float *)malloc(sizeof(float));
		delay->matrix[0][0] = constant_node.attribute("max").as_float();
		parse_delay_input_output(constant_node.attribute("in_port").as_string(), constant_node.attribute("out_port").as_string(), delay);
		print_tabs(level, "%s -> %s Delay: %g\n", constant_node.attribute("in_port").as_string(), constant_node.attribute("out_port").as_string(), delay->matrix[0][0]);

		interconnect->delays.push_back(delay);
	}

/*	for (auto constant_node = interconnect_node.child("delay_constant"); constant_node; constant_node = constant_node.next_sibling("delay_constant")) {*/
/*		Delay *delay = new Delay();*/
/*		delay->type = DELAY_CONSTANT;*/
/*		delay->matrix = (float **)malloc(sizeof(float *));*/
/*		delay->matrix[0] = (float *)malloc(sizeof(float));*/
/*		delay->matrix[0][0] = constant_node.attribute("max").as_float();*/
/*		delay->input = constant_node.attribute("in_port").as_string();*/
/*		delay->output = constant_node.attribute("out_port").as_string();*/
/*		print_tabs(level, "%s -> %s Delay: %g\n", delay->input.c_str(), delay->output.c_str(),delay->matrix[0][0]);*/
/**/
/*		interconnect->delays.push_back(delay);*/
/*	}*/
}

void Arch::parse_interconnect(const xml_node &interconnect_node, Mode *mode, int level)
{
	for (auto direct_node = interconnect_node.child("direct"); direct_node; direct_node = direct_node.next_sibling("direct")) {
		Interconnect *interconnect = new Interconnect();
		auto name = direct_node.attribute("name");
		interconnect->name = name.as_string();	
		print_tabs(level, "Direct: %s\n", interconnect->name.c_str());

		parse_delays(direct_node, interconnect, level+1);

		mode->interconnects.push_back(interconnect);
	}

	for (auto mux_node = interconnect_node.child("mux"); mux_node; mux_node = mux_node.next_sibling("mux")) {
		Interconnect *interconnect = new Interconnect();
		auto name = mux_node.attribute("name");
		interconnect->name = name.as_string();	
		print_tabs(level, "Mux: %s\n", interconnect->name.c_str());

		parse_delays(mux_node, interconnect, level+1);

		mode->interconnects.push_back(interconnect);
	}
	
	for (auto complete_node = interconnect_node.child("complete"); complete_node; complete_node = complete_node.next_sibling("complete")) {
		Interconnect *interconnect = new Interconnect();
		auto name = complete_node.attribute("name");
		interconnect->name = name.as_string();	
		print_tabs(level, "Complete: %s\n", interconnect->name.c_str());

		parse_delays(complete_node, interconnect, level+1);

		mode->interconnects.push_back(interconnect);
	}
}

void Arch::process_lut_class(Pb *lut_pb)
{
	/* Wire mode */
	Mode *wire_mode = new Mode();
	wire_mode->name = "wire";

	Interconnect *wire_mode_interconnect = new Interconnect();

	wire_mode_interconnect->name = "complete:" + lut_pb->name;

	Delay *wire_delay = new Delay();

	wire_delay->type = DELAY_CONSTANT;
	wire_delay->matrix = (float **)malloc(sizeof(float *));
	wire_delay->matrix[0] = (float *)malloc(sizeof(float));
	wire_delay->matrix[0][0] = 235E-12L;

	wire_delay->in_port.instance = lut_pb->name;
	wire_delay->in_port.instance_low = 0;
	wire_delay->in_port.instance_high = INT_MAX;
	wire_delay->in_port.port = "in";
	wire_delay->in_port.port_low = 0;
	wire_delay->in_port.port_high = INT_MAX;

	wire_delay->out_port.instance = lut_pb->name;
	wire_delay->out_port.instance_low = 0;
	wire_delay->out_port.instance_high = INT_MAX;
	wire_delay->out_port.port = "out";
	wire_delay->out_port.port_low = 0;
	wire_delay->out_port.port_high = INT_MAX;

	wire_mode_interconnect->delays.push_back(wire_delay);

	wire_mode->interconnects.push_back(wire_mode_interconnect);

	lut_pb->modes.insert(make_pair(wire_mode->name, wire_mode));

	/* LUT mode */
	Mode *lut_mode = new Mode();
	lut_mode->name = lut_pb->name;

	Interconnect *lut_mode_interconnect_0 = new Interconnect();

	lut_mode_interconnect_0->name = "direct:" + lut_pb->name;

	Delay *lut_delay_0 = new Delay();

	lut_delay_0->type = DELAY_CONSTANT;
	lut_delay_0->matrix = (float **)malloc(sizeof(float *));
	lut_delay_0->matrix[0] = (float *)malloc(sizeof(float));
	lut_delay_0->matrix[0][0] = 0;

	lut_delay_0->in_port.instance = lut_pb->name;
	lut_delay_0->in_port.instance_low = 0;
	lut_delay_0->in_port.instance_high = INT_MAX;
	lut_delay_0->in_port.port = "in";
	lut_delay_0->in_port.port_low = 0;
	lut_delay_0->in_port.port_high = INT_MAX;

	lut_delay_0->out_port.instance = "lut";
	lut_delay_0->out_port.instance_low = 0;
	lut_delay_0->out_port.instance_high = INT_MAX;
	lut_delay_0->out_port.port = "in";
	lut_delay_0->out_port.port_low = 0;
	lut_delay_0->out_port.port_high = INT_MAX;

	lut_mode_interconnect_0->delays.push_back(lut_delay_0);

	Interconnect *lut_mode_interconnect_1 = new Interconnect();

	lut_mode_interconnect_1->name = "direct:" + lut_pb->name;

	Delay *lut_delay_1 = new Delay();

	lut_delay_1->type = DELAY_CONSTANT;
	lut_delay_1->matrix = (float **)malloc(sizeof(float *));
	lut_delay_1->matrix[0] = (float *)malloc(sizeof(float));
	lut_delay_1->matrix[0][0] = 0;

	lut_delay_1->in_port.instance = "lut";
	lut_delay_1->in_port.instance_low = 0;
	lut_delay_1->in_port.instance_high = INT_MAX;
	lut_delay_1->in_port.port = "out";
	lut_delay_1->in_port.port_low = 0;
	lut_delay_1->in_port.port_high = INT_MAX;

	lut_delay_1->out_port.instance = lut_pb->name;
	lut_delay_1->out_port.instance_low = 0;
	lut_delay_1->out_port.instance_high = INT_MAX;
	lut_delay_1->out_port.port = "out";
	lut_delay_1->out_port.port_low = 0;
	lut_delay_1->out_port.port_high = INT_MAX;

	lut_mode_interconnect_1->delays.push_back(lut_delay_1);

	lut_mode->interconnects.push_back(lut_mode_interconnect_0);
	lut_mode->interconnects.push_back(lut_mode_interconnect_1);

	
	Pb *child_lut_pb = new Pb();
	child_lut_pb->name = "lut";
	
	lut_mode->pbs.insert(make_pair(child_lut_pb->name, child_lut_pb));

	lut_pb->modes.insert(make_pair(lut_mode->name, lut_mode));
}

void Arch::parse_pb_type(const xml_node &pb_node, Pb *pb, int level)
{
	auto mode_node = pb_node.child("mode");
	auto pb_name = pb_node.attribute("name");
	int num_child;

	pb->name = pb_name.as_string();
	print_tabs(level, "pb_type name: %s\n", pb->name.c_str());
	pb->blif_model = pb_node.attribute("blif_model").as_string();
	pb->class_name = pb_node.attribute("class").as_string();

	if (pb->blif_model != "") {
		if (pb->class_name == "lut") {
			process_lut_class(pb);
		}
	} else {
		if (mode_node.empty()) {
			Mode *mode = new Mode();
			mode->name = pb_name.as_string();

			print_tabs(level, "Single mode\n");
			num_child = 0;
			for (auto child_node = pb_node.child("pb_type"); child_node; child_node = child_node.next_sibling("pb_type")) {
				Pb *child_pb = new Pb();
				parse_pb_type(child_node, child_pb, level+1);
				mode->pbs.insert(make_pair(child_pb->name, child_pb));
				num_child++;
			}
			print_tabs(level, "Single mode has %d children\n", num_child);

			parse_interconnect(pb_node.child("interconnect"), mode, level);

			pb->modes.insert(make_pair(mode->name, mode));
		} else {
			for (; mode_node; mode_node = mode_node.next_sibling("mode")) {
				auto mode_name = mode_node.attribute("name");
				Mode *mode = new Mode();
				mode->name = mode_name.as_string();
				print_tabs(level, "Mode: %s\n", mode->name.c_str());
				num_child = 0;
				for (auto child_node = mode_node.child("pb_type"); child_node; child_node = child_node.next_sibling("pb_type")) {
					Pb *child_pb = new Pb();
					parse_pb_type(child_node, child_pb, level+1);
					mode->pbs.insert(make_pair(child_pb->name, child_pb));
					num_child++;
				}
				print_tabs(level, "Mode: %s has %d children\n", mode->name.c_str(), num_child);

				parse_interconnect(mode_node.child("interconnect"), mode, level);

				pb->modes.insert(make_pair(mode->name, mode));
			}
			assert(pb_node.child("interconnect").empty());
		}
	}
/*	printf("\n");*/
}

void Arch::parse_complex_block_list(const xml_node &cb_node)
{
	assert(!strcmp(cb_node.name(), "complexblocklist"));

	for (auto pb_node = cb_node.child("pb_type"); pb_node; pb_node = pb_node.next_sibling("pb_type")) {
		Pb *pb = new Pb();
		parse_pb_type(pb_node, pb, 0);
		pbs.insert(make_pair(pb->name, pb));
	}
}

void Arch::load(const string &filename)
{
	xml_document doc;
	xml_parse_result result = doc.load_file(filename.c_str());
	assert(result.status == status_ok);
	xml_node root_node = doc.document_element();

/*	printf("root node name: %s\n", root_node.name());*/
	parse_complex_block_list(root_node.child("complexblocklist"));
}
