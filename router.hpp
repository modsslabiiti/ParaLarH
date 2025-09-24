#ifndef ROUTER_H
#define ROUTER_H

#include "graph.hpp"
#include "basic_types.hpp"

namespace boost {
	enum edge_object_t { edge_object };

	BOOST_INSTALL_PROPERTY(edge, object);
}

typedef struct _GridGraphProperty {
	int nx;
	int ny;
} GridGraphProperty;

typedef struct _GridGraphEdgeObject {
	int util;
	int max_util;
	float mult;
	float delay;
	float slp;
//	int slp;
//	float stp_size;
//	float slp_const;
	//float stp_siz;e_const;
//	float slp_array[1000];
	int m;
	float T1;
	float T2;
//	int util_previous;
//	float mult_previous;
	int violation;
} GridGraphEdgeObject;

typedef property<vertex_name_t, string> GridGraphVertexProperty;

typedef property<edge_object_t, GridGraphEdgeObject,
	property<edge_weight_t, float>> GridGraphEdgeProperty;

typedef Graph<adjacency_list<vecS, vecS, undirectedS, GridGraphVertexProperty, GridGraphEdgeProperty, GridGraphProperty>, Point> GridGraph;

GridGraph create_grid_graph(int nx, int ny);

class Netlist;

class Router {
	public:
		Router(int nx, int ny) : nx(nx), ny(ny), grid(create_grid_graph(nx, ny)) {
		}
		void parallel_route_tbb_new_netlist(Netlist &netlist, int num_threads);
		int nx, ny;
		GridGraph grid;
};

#endif
