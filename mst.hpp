#ifndef MST_H
#define MST_H

#include "graph.hpp"

using namespace std;
using namespace boost;


template<typename CompleteGraph, typename Graph>
struct Triple {
	typename CompleteGraph::vertex_descriptor vertices[3];
	float best_d;
	typename Graph::vertex_descriptor best_steiner;
/*	bool operator<(const Triple &other) const {*/
};

template<typename G>
class Mst {
	public:
		/*		typedef typename graph_traits<G>::vertex_descriptor Vertex;*/
		/*		typedef typename graph_traits<G>::vertex_iterator VertexIterator;*/

		template<typename Container>
			list<typename G::edge_descriptor> kmb(const G &g, const Container &points) const;

		template<typename Container, typename PreviousSteiner>
			list<typename G::edge_descriptor> zel(const G &g, const Container &points, PreviousSteiner &previous_steiners) const;

		template<typename Container, typename PreviousSteiner>
			void parallel_zel(const G &g, const Container &points, PreviousSteiner &previous_steiners) const;
		template<typename Container, typename PreviousSteiner>
			void parallel_zel_new(const G &g, const Container &points, PreviousSteiner &previous_steiners) const;
		template<typename Container, typename PreviousSteiner>
			void parallel_zel_new_reduce(const G &g, const Container &points, PreviousSteiner &previous_steiners) const;
/*		template<typename FilteredGraph, typename PredMap, typename SteinerVertices>*/
/*			void remove_non_steiner_points(const FilteredGraph &g, PredMap &pred_map, const SteinerVertices &steiner_vertices);*/

		vector<vector<float>> all_dist;
		vector<vector<typename G::vertex_descriptor>> all_pred;
		
		void init(const G &g);
};

template<class G, class PredMap>
struct mst_edge_filter {
	mst_edge_filter() : g(NULL), pred_map(NULL) {}
	mst_edge_filter(const G *g, const PredMap *pred_map) : g(g), pred_map(pred_map) {}

	template<class Edge>
	bool operator()(const Edge &e) const {
/*		printf("src: %d tar: %d pred: %d\n", source(e, *g), target(e, *g), (*pred_map)[target(e, *g)]);*/
		return (*pred_map)[target(e, *g)] == source(e, *g) || (*pred_map)[source(e, *g)] == target(e, *g);
	}

	const G *g;
	const PredMap *pred_map;
};

template<class G, class PredMap>
struct mst_vertex_filter {
	mst_vertex_filter() : g(NULL), pred_map(NULL) {}
	mst_vertex_filter(const G *g, const PredMap *pred_map) : g(g), pred_map(pred_map) {}

	template<class Vertex>
	bool operator()(const Vertex &v) const {
		return (*pred_map)[v] != v;
	}

	const G *g;
	const PredMap *pred_map;
};

#include "mst_impl.hpp"

#endif
