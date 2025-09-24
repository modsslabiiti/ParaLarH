#ifndef TIMING_H
#define TIMING_H

namespace boost {
	enum vertex_t_arr_t { vertex_t_arr };

	BOOST_INSTALL_PROPERTY(vertex, t_arr);
}

template<typename RouteTree, typename ColorMap>
void get_delay(RouteTree &route_tree, typename RouteTree::vertex_descriptor src, float delay, ColorMap &visited)
{
	typename RouteTree::out_edge_iterator e, e_end;

	assert(src != graph_traits<RouteTree>::null_vertex());

	put(vertex_t_arr, route_tree, src, delay);
/*	printf("level: %d\n", delay);*/
	visited[src] = true;

	for (std::tie(e, e_end) = out_edges(src, route_tree); e != e_end; ++e) {
		assert(source(*e, route_tree) == src);
		if (!visited[target(*e, route_tree)]) {
			get_delay(route_tree, target(*e, route_tree), delay + get(edge_weight, route_tree, *e), visited);
		}
	}
}


template<typename TimingGraph>
float get_critical_path_delay(TimingGraph &tgraph, float *average_delay, struct timing_node *crit_tnode)
{
/*	assert(is_connected(tgraph));*/
	adjacency_list<vecS, vecS, undirectedS> undirected_tgraph;
	for (int i = 0; i < num_vertices(tgraph); ++i) {
		assert(add_vertex(undirected_tgraph) < num_vertices(tgraph));
	}
	typename TimingGraph::edge_iterator tei, tei_end;
	for (std::tie(tei, tei_end) = edges(tgraph); tei != tei_end; ++tei) {
		add_edge(source(*tei, tgraph), target(*tei, tgraph), undirected_tgraph);
	}

	vector<int> comp(num_vertices(undirected_tgraph));
	auto comp_map = make_iterator_property_map(comp.begin(), get(vertex_index, undirected_tgraph));
/*	printf("number of connected components: %d\n", connected_components(undirected_tgraph, comp_map));*/
/*	bool connected = true;*/
	for (int i = 0; i < comp.size(); ++i) {
/*		if (comp[i] != 0) {*/
/*			printf("[%d] Block: %s %s Port: %s Pin: %d\n", comp[i], tgraph[i].pin->port->block->name.c_str(), tgraph[i].pin->port->block->instance.c_str(), tgraph[i].pin->port->name.c_str(), tgraph[i].pin->index);*/
/*		}*/
	}

	vector<typename TimingGraph::vertex_descriptor> res;
	topological_sort(tgraph, back_inserter(res));

	float critical_path_delay = FLT_MIN;
	float total_delay = 0;
	int count = 0;

	typename TimingGraph::vertex_iterator vi, vi_end;

	for (std::tie(vi, vi_end) = vertices(tgraph); vi != vi_end; ++vi) {
		tgraph[*vi].arrival_time = FLT_MIN;
		tgraph[*vi].required_time = FLT_MAX;
		tgraph[*vi].prev = nullptr;
	}

	/* forward propagation */
	for (auto iter = res.rbegin(); iter != res.rend(); ++iter) {
		if (in_degree(*iter, tgraph) == 0) {
			tgraph[*iter].arrival_time = 0;
/*			printf("Input: %s %s %s %d Delay: %f\n", tgraph[*iter].pin->port->block->name.c_str(), tgraph[*iter].pin->port->block->instance.c_str(), tgraph[*iter].pin->port->name.c_str(), tgraph[*iter].pin->index, tgraph[*iter].arrival_time);*/
		}

		/* track primary output arrival time */
		if (out_degree(*iter, tgraph) == 0) {
/*			critical_path_delay = std::max(critical_path_delay, tgraph[*iter].arrival_time);*/
			if (tgraph[*iter].arrival_time > critical_path_delay) {
				critical_path_delay = tgraph[*iter].arrival_time;
				*crit_tnode = tgraph[*iter];
			}
			if (tgraph[*iter].arrival_time > 0) {
				total_delay += tgraph[*iter].arrival_time;
				count++;
			}
/*			if (tgraph[*iter].arrival_time > 0) {*/
/*			printf("Output: %s %s %s %d Delay: %f\n", tgraph[*iter].pin->port->block->name.c_str(), tgraph[*iter].pin->port->block->instance.c_str(), tgraph[*iter].pin->port->name.c_str(), tgraph[*iter].pin->index, tgraph[*iter].arrival_time);*/
/*			}*/
		}

		typename TimingGraph::out_edge_iterator ei, ei_end;

		for (std::tie(ei, ei_end) = out_edges(*iter, tgraph); ei != ei_end; ++ei) {
			assert(source(*ei, tgraph) == *iter);
			if (tgraph[*iter].arrival_time + get(edge_weight, tgraph, *ei) > tgraph[target(*ei, tgraph)].arrival_time) {
				tgraph[target(*ei, tgraph)].arrival_time = tgraph[*iter].arrival_time + get(edge_weight, tgraph, *ei);
				tgraph[target(*ei, tgraph)].prev = &tgraph[*iter];
			}

/*			tgraph[target(*ei, tgraph)].arrival_time = std::max(, );*/
		}
	}

	/* backward propagation */
	for (auto iter = res.begin(); iter != res.end(); ++iter) {
		if (out_degree(*iter, tgraph) == 0) {
			tgraph[*iter].required_time = critical_path_delay;
		}

		typename TimingGraph::in_edge_iterator ei, ei_end;
		 
		for (std::tie(ei, ei_end) = in_edges(*iter, tgraph); ei != ei_end; ++ei) {
			tgraph[source(*ei, tgraph)].required_time = std::min(tgraph[source(*ei, tgraph)].required_time, tgraph[*iter].required_time - get(edge_weight, tgraph, *ei));
		}
	}

	*average_delay = total_delay/count;

	return critical_path_delay;

}
#endif
