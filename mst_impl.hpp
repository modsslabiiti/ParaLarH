#ifndef MST_IMPL_H
#define MST_IMPL_H

#include "mst.hpp"
#include "router.hpp"

/*#define VERBOSE*/
/*#define PRINT_MST_BEFORE*/
/*#define PRINT_MST_AFTER*/
/*#define PRINT_STATS*/

template<typename G>
void Mst<G>::init(const G &g) 
{
	all_dist.resize(num_vertices(g), vector<float>(num_vertices(g)));
	all_pred.resize(num_vertices(g), vector<typename G::vertex_descriptor>(num_vertices(g)));

/*	printf("num_vertices: %d\n", num_vertices(g));*/

	for (const auto &v : g.vertices()) {
		int v_index = get(vertex_index, g, v);
/*		printf("v_index: %d, v: %d\n", v_index, v);*/
		assert(v_index < all_dist.size() && v == v_index);
		dijkstra_shortest_paths(g, v, 
				weight_map(get(edge_weight, g)).
				distance_map(make_iterator_property_map(all_dist[v_index].begin(), get(vertex_index, g))).
				predecessor_map(make_iterator_property_map(all_pred[v_index].begin(), get(vertex_index, g))));
	}
}

/*template<typename G>*/
/*template<typename FilteredGraph, typename PredMap, typename SteinerVertices>*/
/*void Mst<G>::remove_non_steiner_points(const FilteredGraph &g, PredMap &pred_map, const SteinerVertices &steiner_vertices)*/
/*{*/
/*	typename FilteredGraph::vertex_iterator v, v_end;*/
/*	for (tie(v, v_end) = vertices(g); v != v_end; ++v) {*/
/*		typename FilteredGraph::vertex_descriptor cur_v = *v;*/
/*		while (out_degree(cur_v, g) == 1 && steiner_vertices.find(cur_v) == steiner_vertices.end()) {*/
/*			typename FilteredGraph::adjacency_iterator n, n_end;*/
/*			tie(n, n_end) = adjacent_vertices(cur_v, g);*/
/*			assert(next(n, 1) == n_end);*/
/*			pred_map[cur_v] = cur_v;*/
/*			cur_v = *n;*/
/*		}*/
/*	}*/
/*}*/

template<class Graph>
struct vertex_label_writer {
	const Graph *g;
	vertex_label_writer(const Graph *g) : g(g) {}
	template<class VertexOrEdge>
	void operator()(std::ostream& out, const VertexOrEdge& v) const {
		out << "[label=" << escape_dot_string(get(vertex_name, *g, v)) << "]";
	}
};

template<typename G>
template<typename Container>
list<typename G::edge_descriptor> Mst<G>::kmb(const G &g, const Container &points) const
{
	typedef Graph<adjacency_list<vecS, vecS, undirectedS, property<vertex_name_t, string>, property<edge_weight_t, float>>, typename G::object_type> TempGraph;
	TempGraph dist_g;
#ifdef VERBOSE
	printf("kmb input: ");
	for (const auto &p : points) {
		printf("%d,%d ", p.x, p.y);
/*		dist_g.add_vertex(p);*/
	}
	printf("\n");
#endif

/*	if (do_init) {*/
/*		init(g);*/
/*	}*/

	for (auto p1 = points.begin(); p1 != points.end(); ++p1) {
		for (auto p2 = std::next(p1, 1); p2 != points.end(); ++p2) {
			typename TempGraph::edge_descriptor e = dist_g.add_edge(*p1, *p2);
			put(edge_weight, dist_g, e, all_dist[g.vertex(*p1)][g.vertex(*p2)]);
#ifdef VERBOSE
			cout << "[KMB] " <<  *p1 << " -> " << *p2 << " weight: " << all_dist[g.vertex(*p1)][g.vertex(*p2)] << endl;
#endif
		}
	}
	{
		int n = num_vertices(dist_g);
/*		if (n != points.size()) {*/
/*			printf("n: %d points: %d\n", n, points.size());*/
/*		}*/
		assert(n == points.size());
/*		printf("[KMB] g_metric num_edges: %d\n", n);*/
		assert(n > 0);
		assert(num_edges(dist_g) == n*(n-1)/2);
	}
	
	vector<typename TempGraph::vertex_descriptor> mst_pred(num_vertices(dist_g));
	auto mst_pred_map = make_iterator_property_map(mst_pred.begin(), get(vertex_index, dist_g));

	prim_minimum_spanning_tree(dist_g, mst_pred_map, weight_map(get(edge_weight, dist_g)));

	TempGraph shortest_path_g;

	for (const auto &v : dist_g.vertices()) {
		typename TempGraph::vertex_descriptor src = mst_pred_map[v];
		typename TempGraph::vertex_descriptor dst = v;
		typename G::vertex_descriptor src_g = g.vertex(dist_g.object(src));
		typename G::vertex_descriptor dst_g = g.vertex(dist_g.object(dst));
		if (src != dst) {
#ifdef VERBOSE
			cout << "[KMB] Replacing " << dist_g.object(src) << " -> " << dist_g.object(dst) << "with shortest path\n";
#endif
			const auto &pred = all_pred[src_g];
			auto p = dst_g;
#ifdef VERBOSE
			printf("(%d,%d)\n", g.object(p).x, g.object(p).y);
#endif
			while (p != src_g) {
/*				assert(abs(g.object(pred[p]).x-g.object(p).x) == 1 || abs(g.object(pred[p]).y-g.object(p).y) == 1);*/
/*				g.edge(pred[p], p)*/
				typename TempGraph::edge_descriptor e = shortest_path_g.add_edge(g.object(pred[p]), g.object(p));

				put(edge_weight, shortest_path_g, e, all_dist[pred[p]][p]);

				put(vertex_name, shortest_path_g, shortest_path_g.vertex(g.object(pred[p])), get(vertex_name, g, pred[p]));
				put(vertex_name, shortest_path_g, shortest_path_g.vertex(g.object(p)), get(vertex_name, g, p));

				p = pred[p];
#ifdef VERBOSE
				printf("(%d,%d)\n", g.object(p).x, g.object(p).y);
#endif
			}
		}
	}

	extern string current_net_details;
/*	ofstream shortest_path_g_file("output/" + current_net_details + "kmbspg");*/
/*	write_graphviz(shortest_path_g_file, shortest_path_g, vertex_label_writer<TempGraph>(&shortest_path_g));*/

/*	assert(num_vertices(shortest_path_g)-1 == num_edges(shortest_path_g));*/

	mst_pred.resize(num_vertices(shortest_path_g));
	mst_pred_map = make_iterator_property_map(mst_pred.begin(), get(vertex_index, shortest_path_g));
	prim_minimum_spanning_tree(shortest_path_g, mst_pred_map, weight_map(get(edge_weight, shortest_path_g)));
/*	for (const auto &v : shortest_path_g.vertices()) {*/
/*		if (mst_pred_map[v] != v) {*/
/*			cout << "Final edge: " << shortest_path_g.object(mst_pred_map[v]) << " -> " << shortest_path_g.object(v) << endl;*/
/*		}*/
/*	}*/

/*	typedef filtered_graph<TempGraph, mst_edge_filter<TempGraph, decltype(mst_pred_map)>> FilteredGraph;*/
/*	FilteredGraph shortest_path_g_mst(shortest_path_g, mst_edge_filter<TempGraph, PredMap>(&shortest_path_g, &mst_pred_map));*/
/*	typename FilteredGraph::edge_iterator e, e_end;*/
/*	for (tie(e, e_end) = edges(shortest_path_g_mst); e != e_end; ++e) {*/
/*		cout << "Final edge: " << shortest_path_g.object(source(*e, shortest_path_g)) << " -> " << shortest_path_g.object(target(*e, shortest_path_g)) << endl;*/
/*	}*/
	list<typename TempGraph::edge_descriptor> edges_to_remove;
	for (const auto &e : shortest_path_g.edges()) {
/*		cout << "Final edge: " << shortest_path_g.object(source(e, shortest_path_g)) << " -> " << shortest_path_g.object(target(e, shortest_path_g)) << endl;*/
/*		cout << "mst_pred[source] " << shortest_path_g.object(mst_pred_map[source(e, shortest_path_g)]) << " mst_pred[target] " << shortest_path_g.object(mst_pred_map[target(e, shortest_path_g)]) << endl;*/
		if (mst_pred_map[source(e, shortest_path_g)] != target(e, shortest_path_g) && mst_pred_map[target(e, shortest_path_g)] != source(e, shortest_path_g)) {
/*			cout << "Removing\n";*/
			edges_to_remove.push_back(e);	
		}
	}
	for (const auto &e : edges_to_remove) {
		remove_edge(e, shortest_path_g);
	}

/*	ofstream shortest_path_g_mst_file("output/" + current_net_details + "kmbspgmst");*/
/*	write_graphviz(shortest_path_g_mst_file, shortest_path_g, vertex_label_writer<TempGraph>(&shortest_path_g));*/

	assert(is_tree(shortest_path_g));

	/* fast lookup to check whether a point is steiner */
	set<typename TempGraph::vertex_descriptor> steiner_vertices;
	for (const auto &p : points) {
		steiner_vertices.insert(shortest_path_g.vertex(p));
	}
	TempGraph shortest_path_mst;
	for (const auto &v : shortest_path_g.vertices()) {
		/* if non-steiner leaf */
		typename TempGraph::vertex_descriptor cur_v = v;
		while (out_degree(cur_v, g) == 1 && steiner_vertices.find(cur_v) == steiner_vertices.end()) {
			typename TempGraph::adjacency_iterator n, n_end;
			std::tie(n, n_end) = adjacent_vertices(cur_v, shortest_path_g);
			assert(next(n, 1) == n_end);
			assert(mst_pred_map[cur_v] == *n || mst_pred_map[*n] == v);
			remove_edge(v, *n, shortest_path_g);
			/* remove the edge */

/*				if (mst_pred_map[cur_v] == *n) {*/
/*					mst_pred_map[cur_v] = cur_v;*/
/*				} else {*/
/*					mst_pred_map[*n] = *n;*/
/*				}*/
			cur_v = *n;
		}
	}
/*	remove_non_steiner_points(shortest_path_g_mst, mst_pred_map, points_lookup);*/
	assert(is_tree(shortest_path_g));

	list<typename G::edge_descriptor> output;
	for (const auto &e : shortest_path_g.edges()) {
/*		if (mst_pred_map[v] != v) {*/
			output.push_back(g.edge(shortest_path_g.object(source(e, shortest_path_g)), shortest_path_g.object(target(e, shortest_path_g))));
/*		}*/
	}

	typedef Graph<adjacency_list<vecS, vecS, undirectedS, property<vertex_name_t, string>>, typename G::object_type> VerifyGraph;
	VerifyGraph verify;
	for (const auto &e : output) {
		auto u = verify.add_vertex(g.object(source(e, g)));
		auto v = verify.add_vertex(g.object(target(e, g)));
		put(vertex_name, verify, u, get(vertex_name, g, source(e, g)));
		put(vertex_name, verify, v, get(vertex_name, g, target(e, g)));
		verify.add_edge(u, v);
	}
	/*verification*/
	assert(is_tree(verify));
	for (const auto &p : points) {
		auto v = verify.vertex(p);
		assert(v >= 0 && v < num_vertices(verify));
	}

	return output;
}

template<typename Tri>
struct triple_comp {
	bool operator()(const Tri &t1, const Tri &t2) const {
		return t1->best_d > t2->best_d;
	}
};

template<typename G>
struct edge_weight_filter {
	typedef typename G::edge_descriptor edge;
	typedef typename property_map<G, edge_weight_t>::const_type edge_property_map;
	typedef typename property_traits<edge_property_map>::value_type value_type;
	typedef typename property_traits<edge_property_map>::key_type key_type;
	typedef typename property_traits<edge_property_map>::category category;
	typedef typename property_traits<edge_property_map>::reference reference;

	const edge_property_map &prop;
	edge uv;
	edge vw;

	edge_weight_filter(const edge_property_map &prop, edge uv, edge vw) : prop(prop), uv(uv), vw(vw) {
	}
};

template<typename G>
typename property_traits<typename edge_weight_filter<G>::edge_property_map>::value_type get(const edge_weight_filter<G> &filter, typename edge_weight_filter<G>::edge e)
{
/*	return (e == filter.uv || e == filter.vw) ? 0 : get(filter.prop, e);*/
	if (e == filter.uv || e == filter.vw) {
/*		printf("returning 0\n");*/
		return 0;
	} else { 
		float lol = get(filter.prop, e);
/*		printf("returning %f\n", lol);*/
		return lol;
	}
}

template<typename G>
template<typename Container, typename PreviousSteiner>
void Mst<G>::parallel_zel_new_reduce(const G &g, const Container &input_points, PreviousSteiner &previous_steiners) const
{
/*	if (do_init) {*/
/*		init(g);*/
/*	}*/

/*	VertexIterator v1, v1_end;*/
/*	VertexIterator v2, v2_end;*/

/*	for (tie(v1, v1_end) = vertices(g); v != v1_end; ++v1) {*/
/*		for (tie(v2, v2_end) = vertices(g); v2 != v2_end; ++v2) {*/
/*			Point p1 = g.object(*v);*/
/*			Point p2 = g.object(*v2);*/
/*			assert(all_dist[*v][*v2] == abs(p1.x - p2.x) + abs(p1.y - p2.y));*/
/*		}*/
/*	}*/

	typedef Graph<adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, float>>, typename G::object_type> TempGraph;
	TempGraph f;

#ifdef VERBOSE
	printf("zel input: ");
	for (const auto &p : points) {
		printf("%d,%d ", p.x, p.y);
/*		f.add_vertex(p);*/
	}
	printf("\n");
#endif
	set<typename G::object_type> point_set;
	
	for (const auto &p : input_points) {
		point_set.insert(p);
	}
	vector<typename G::object_type> points;
	for (const auto &p : point_set) {
		points.emplace_back(p);
	}

	/* construct complete graph */
	for (auto p1 = points.begin(); p1 != points.end(); ++p1) {
		for (auto p2 = std::next(p1, 1); p2 != points.end(); ++p2) {
			typename TempGraph::edge_descriptor e;
			e = f.add_edge(*p1, *p2);
			put(edge_weight, f, e, all_dist[g.vertex(*p1)][g.vertex(*p2)]);
		}
	}
	int n = points.size();
	assert(num_vertices(f) == n && num_edges(f) == n*(n-1)/2);

	/* create all possible triples */
	vector<std::unique_ptr<Triple<TempGraph, G>>> triples;
	int num_combinations = std::max(n*(n-1)*(n-2)/6, 0);
	triples.reserve(num_combinations);
	int comb = 0;

	boost::timer::cpu_timer timer;
	timer.start();

	for (int i = 0; i < n; ++i) {
		for (int j = i+1; j < n; ++j) {
			for (int k = j+1; k < n; ++k) {
				Triple<TempGraph, G> *t = new Triple<TempGraph, G>();
				t->vertices[0] = f.vertex(points[i]);
				t->vertices[1] = f.vertex(points[j]);
				t->vertices[2] = f.vertex(points[k]);

#ifdef VERBOSE
				printf("i[%d]=(%d,%d)[%lu], k[%d]=(%d,%d)[%lu], j[%d]=(%d,%d)[%lu]\n", 
						i, points[i].x, points[i].y, t->vertices[0],
						j, points[j].x, points[j].y, t->vertices[1],
						k, points[k].x, points[k].y, t->vertices[2]);
#endif

				triples.emplace_back(t);
/*				push_heap(triples.begin(), triples.end(), triple_comp<Triple<TempGraph, G> *>());*/
				comb++;
			}
		}
	}
	assert(comb == num_combinations);

/*	tbb::affinity_partitioner ap;*/

	/* previous grain size: 128 */
	tbb::parallel_for(tbb::blocked_range<int>(0, triples.size()), [&g,&f,&points,&triples,this](const tbb::blocked_range<int> &range) -> void {
		for (int i = range.begin(); i != range.end(); ++i) {
			float min_d = FLT_MAX;
			typename G::vertex_descriptor min_v = G::null_vertex();

			for (const auto &v : g.vertices()) {
				float total_d = 0;
				for (int j = 0; j < 3; ++j) {
					total_d += all_dist[g.vertex(f.object(triples[i]->vertices[j]))][v];
				}
#ifdef VERBOSE
				printf("min: %f current: %f\n", min_d, total_d);
				for (int lol = 0; lol < 3; ++lol) {
					printf("indices: %d ", indices[lol]);
				}
				printf("\n");
#endif
				if (total_d < min_d) {
					min_d = total_d;
					min_v = v;
				}
			}

			assert(min_v != num_vertices(g));
			triples[i]->best_d = min_d;
			triples[i]->best_steiner = min_v;
#ifdef VERBOSE
			printf("best: %f (%d,%d)\n", t->best_d, t->best_steiner.x, t->best_steiner.y);
#endif
		}
	});

	timer.stop();
	boost::timer::cpu_times triple_gen_elapsed = timer.elapsed();
	boost::timer::nanosecond_type triple_gen_time = triple_gen_elapsed.system + triple_gen_elapsed.user;

	/* main algo loop */
	bool done = false;
	set<typename G::object_type> steiners;
	vector<typename TempGraph::vertex_descriptor> mst_pred(num_vertices(f));
	auto mst_pred_map = make_iterator_property_map(mst_pred.begin(), get(vertex_index, f));
	map<typename TempGraph::edge_descriptor, float> contracted_edges;

	timer.start();
	while (!done) {

		prim_minimum_spanning_tree(f, mst_pred_map, weight_map(get(edge_weight, f)));
/*		mst_edge_filter<PredMap> filter(mst_pred_map);*/
/*		filtered_graph<TempGraph, mst_edge_filter<PredMap>> */
		float mst_f_total_before = 0;
		map<typename TempGraph::edge_descriptor, float> mst_before_edges;
#ifdef PRINT_MST_BEFORE 
		printf("MST Before\n");
#endif
		for (const auto &v : f.vertices()) {
			typename TempGraph::vertex_descriptor pre = mst_pred_map[v];
			if (pre != v) {
#ifdef PRINT_MST_BEFORE
			printf("MST edge (%d,%d) -> (%d, %d)\n", f.object(pre).x, f.object(pre).y, f.object(v).x, f.object(v).y);
#endif
				float w = get(edge_weight, f, f.edge(pre, v));
				mst_f_total_before += w;
				assert(mst_before_edges.find(f.edge(pre, v)) == mst_before_edges.end());
				mst_before_edges[f.edge(pre, v)] = w;
			}
		}

		int pruned = 0;

/*		for (const auto &t : triples) {*/
		/* previous grain size: 1024 */
		auto result = tbb::parallel_deterministic_reduce(tbb::blocked_range<int>(0, triples.size()), std::pair<float, const Triple<TempGraph, G> *>(0, nullptr),
			[&f,&mst_f_total_before,&triples](const tbb::blocked_range<int> &range, const std::pair<float, const Triple<TempGraph, G> *> &initial) -> std::pair<float, const Triple<TempGraph, G> *> {
/*			printf("f:%X mst_f_total_before:%X triples:%X best_triple:%X max_win:%X g_mutex:%X\n", &f, &mst_f_total_before, &triples, &best_triple, &max_win, &g_mutex);*/
/*			TempGraph *f_local_ptr = new TempGraph(f);*/
/*			TempGraph &f_local = *f_local_ptr;*/
			const TempGraph &f_local = f;
/*			const TempGraph f_local(f);*/
/*			vector<typename TempGraph::vertex_descriptor> *mst_pred_local_ptr = new vector<typename TempGraph::vertex_descriptor>(num_vertices(f_local));*/
/*			vector<typename TempGraph::vertex_descriptor> &mst_pred_local = *mst_pred_local_ptr;*/
			vector<typename TempGraph::vertex_descriptor> mst_pred_local(num_vertices(f_local));
			auto mst_pred_map_local = make_iterator_property_map(mst_pred_local.begin(), get(vertex_index, f_local));
			
			float local_max_win = initial.first;
			const Triple<TempGraph, G> *local_best_triple = initial.second;

/*			printf("range begin: %d range end: %d\n", range.begin(), range.end());*/

			for (int i = range.begin(); i != range.end(); ++i) {

				const Triple<TempGraph, G> *t = triples[i].get();

				typename TempGraph::vertex_descriptor u = t->vertices[0];	
				typename TempGraph::vertex_descriptor v = t->vertices[1];	
				typename TempGraph::vertex_descriptor w = t->vertices[2];	

				typename TempGraph::edge_descriptor uv = f_local.edge(u, v);
				typename TempGraph::edge_descriptor vw = f_local.edge(v, w);

/*				float old_uv_length = get(edge_weight, f_local, uv);*/
/*				float old_vw_length = get(edge_weight, f_local, vw);*/

				/* A serious bug in the previous commit caused a large amount of triples to be pruned */
				/* Bug: after pruning a triple, the heap is not updated. therefore, front always remain the same */

#ifdef VERBOSE
				printf("Contracting (%d,%d)[%lu] -> (%d,%d)[%lu] -> (%d,%d)[%lu] old_uv: %f old_vw: %f\n", 
						f_local.object(u).x, f_local.object(u).y, u,
						f_local.object(v).x, f_local.object(v).y, v,
						f_local.object(w).x, f_local.object(w).y, w,
						old_uv_length, old_vw_length);
#endif
				auto filter = edge_weight_filter<TempGraph>(get(edge_weight, f_local), uv, vw);
				prim_minimum_spanning_tree(f_local, mst_pred_map_local, weight_map(filter));

#ifdef PRINT_MST_AFTER 
				printf("MST After\n");
#endif
				typename TempGraph::vertex_iterator vert;
				float mst_f_total_after = 0;
				for (const auto &v : f_local.vertices()) {
					typename TempGraph::vertex_descriptor pre = mst_pred_map_local[v];
					if (pre != v) {
#ifdef PRINT_MST_AFTER 
						printf("MST edge (%d,%d) -> (%d, %d)\n", f_local.object(pre).x, f_local.object(pre).y, f_local.object(v).x, f_local.object(v).y);
#endif
						mst_f_total_after += get(filter, f_local.edge(pre, v));
					}
				}

				float win = mst_f_total_before - mst_f_total_after - t->best_d;

	/*				g_mutex.lock();*/
				if (win > local_max_win) {
					local_best_triple = t;
					local_max_win = win;
				}
/*				g_mutex.unlock();*/
				/*				float offset = 0;*/
				/*				auto c_uv = contracted_edges.find(uv);*/
				/*				if (c_uv != contracted_edges.end()) {*/
				/*					offset += -c_uv->second;*/
				/*				}*/
				/*				auto c_vw = contracted_edges.find(vw);*/
				/*				if (c_vw != contracted_edges.end()) {*/
				/*					offset += -c_vw->second;*/
				/*				}*/
				/*			map<TempGraph::edge_descriptor, float>::iterator b_uv = mst_before_edges.find(uv);*/
				/*			if (b_uv != mst_before_edges.end()) {*/
				/*				offset = b_uv->second;*/
				/*			}*/
				/*			map<TempGraph::edge_descriptor, float>::iterator b_vw = mst_before_edges.find(vw);*/
				/*			if (b_vw != mst_before_edges.end()) {*/
				/*				offset = b_vw->second;*/
				/*			}*/

#ifdef PRINT_STATS
				printf("mst_before: %f mst_after: %f offset: %f old_uv: %f old_vw: %f\n", mst_f_total_before, mst_f_total_after, offset, old_uv_length, old_vw_length);
#endif

				/*always true*/
				/*				if (mst_f_total_before - mst_f_total_after > old_uv_length+old_vw_length + offset) {*/
				/*					printf("mst_before: %f mst_after: %f old_uv: %f old_vw: %f offset: %f c_uv: %f c_vw: %f\n", mst_f_total_before, mst_f_total_after, old_uv_length, old_vw_length, offset, c_uv->second, c_vw->second);*/
				/*				}*/
				/*				assert(mst_f_total_before - mst_f_total_after <= old_uv_length+old_vw_length + offset); */

				/* false */
				/*			assert(mst_f_total_before-mst_f_total_after >= old_uv_length || mst_f_total_before-mst_f_total_after >= old_vw_length);*/
				/*			assert(mst_f_total_before - mst_f_total_after + offset == old_uv_length+old_vw_length);*/
				/*			assert(mst_f_total_before - mst_f_total_after == old_uv_length+old_vw_length || */
				/*					mst_f_total_before - mst_f_total_after == old_uv_length ||*/
				/*					mst_f_total_before - mst_f_total_after == old_vw_length );*/
#ifdef VERBOSE
				printf("Total before: %f Total after: %f Best d: %f ", mst_f_total_before, mst_f_total_after, t->best_d);
				printf("Win: %f Max win: %f ", win, max_win);
				if (old_uv_length + old_vw_length < t->best_d) {
					pruned++;
					if (win > 0) {
						printf("WEIRD! win=%f when bound is satisfied\n", win);
					}
					/*					assert(win <= 0);*/
					printf("PRUNED");
				}
				printf("\n");

				if (mst_f_total_before - mst_f_total_after > old_uv_length+old_vw_length) {
					printf("BOUND VIOLATION=%f: mst_before: %f mst_after: %f old_uv: %f old_vw: %f\n", mst_f_total_before-mst_f_total_after-(old_uv_length+old_vw_length), mst_f_total_before, mst_f_total_after, old_uv_length, old_vw_length);
				}
				/*				assert(!(mst_f_total_before - mst_f_total_after > old_uv_length + old_vw_length)); */

#endif
				/*			printf("Win: %f best_d: %f\n", win, t->best_d);*/
			}
/*			delete mst_pred_local_ptr;*/
/*			delete f_local_ptr;*/
			return std::make_pair(local_max_win, local_best_triple);
		},
		[](const std::pair<float, const Triple<TempGraph, G> *> x, const std::pair<float, const Triple<TempGraph, G> *> y) -> std::pair<float, const Triple<TempGraph, G> *> {
			if (x.first > y.first) {
				return x;
			} else {
				return y;
			}
		}
		);

#ifdef VERBOSE
		printf("Pruned: %d num_triples: %d\n", pruned, triples.size());
#endif
		const Triple<TempGraph, G> *best_triple = result.second;
		if (!best_triple) {
			done = true;
		} else {
			typename TempGraph::vertex_descriptor u = best_triple->vertices[0];
			typename TempGraph::vertex_descriptor v = best_triple->vertices[1];
			typename TempGraph::vertex_descriptor w = best_triple->vertices[2];

			typename TempGraph::edge_descriptor uv = f.edge(u, v);
			typename TempGraph::edge_descriptor vw = f.edge(v, w);

			assert(contracted_edges.find(uv) == contracted_edges.end() && contracted_edges.find(vw) == contracted_edges.end());

			contracted_edges[uv] = get(edge_weight, f, uv);
			contracted_edges[vw] = get(edge_weight, f, uv);

			put(edge_weight, f, uv, 0);
			put(edge_weight, f, vw, 0);


/*			printf("Best steiner: (%d,%d)\n", best_triple->best_steiner.x, best_triple->best_steiner.y);*/

			steiners.insert(g.object(best_triple->best_steiner));
		}	
	}
	
	timer.stop();
	boost::timer::cpu_times main_elapsed = timer.elapsed();
	boost::timer::nanosecond_type main_time = main_elapsed.system + main_elapsed.user;

/*	printf("ended\n");*/

#ifdef VERBOSE
	printf("Steiners[count=%d]: ", steiners.size());
	for (const auto &s : steiners) {
		printf("%d,%d ", s.x, s.y);
	}
	printf("\n");
#endif

	if (!previous_steiners.empty()) {
#ifdef VERBOSE
		printf("Previous steiners[count=%d]: ", previous_steiners.size());
		for (const auto &prev_s : previous_steiners) {
			printf("%d,%d ", prev_s.x, prev_s.y);
		}
		printf("\n");
#endif
		for (const auto &s : steiners) {
			float min_d = FLT_MAX;
			const typename G::object_type *min_v = nullptr;
			for (const auto &prev_s : previous_steiners) {
				float d = all_dist[g.vertex(s)][g.vertex(prev_s)];
				if (d < min_d) {
					min_d = d;
					min_v = &prev_s;
				}
			}
#ifdef VERBOSE
			cout << s << " is closest to " << *min_v << " dist=" << min_d << endl;
#endif
		}
		//previous_steiners.clear();
	}
/*	for (const auto &s : steiners) {*/
/*		previous_steiners.push_back(s);*/
/*	}*/
/*	printf("done\n");*/

	std::copy(steiners.begin(), steiners.end(), std::inserter(previous_steiners, previous_steiners.begin()));

/*	steiners.insert(points.begin(), points.end());*/

/*	timer.start();*/
/*	list<typename G::edge_descriptor> tree_edges = kmb(g, steiners);*/
/*	timer.stop();*/
/**/
/*	boost::timer::cpu_times kmb_elapsed = timer.elapsed();*/
/*	boost::timer::nanosecond_type kmb_time = kmb_elapsed.system + kmb_elapsed.user;*/
/**/
/*#ifdef VERBOSE*/
/*	printf("Triple gen: %fms Main: %fms Kmb: %fms\n", (double)triple_gen_time/1000000, (double)main_time/1000000, (double)kmb_time/1000000);*/
/*#endif*/

/*	if (!connected) {*/
/*		ofstream file("bug.viz");*/
/*		write_graphviz(file, verify, vertex_label_writer<VerifyGraph>(&verify));*/
/*		exit(0);*/
/*	}*/
}

template<typename G>
template<typename Container, typename PreviousSteiner>
void Mst<G>::parallel_zel_new(const G &g, const Container &points, PreviousSteiner &previous_steiners) const
{
/*	if (do_init) {*/
/*		init(g);*/
/*	}*/

/*	VertexIterator v1, v1_end;*/
/*	VertexIterator v2, v2_end;*/

/*	for (tie(v1, v1_end) = vertices(g); v != v1_end; ++v1) {*/
/*		for (tie(v2, v2_end) = vertices(g); v2 != v2_end; ++v2) {*/
/*			Point p1 = g.object(*v);*/
/*			Point p2 = g.object(*v2);*/
/*			assert(all_dist[*v][*v2] == abs(p1.x - p2.x) + abs(p1.y - p2.y));*/
/*		}*/
/*	}*/

	typedef Graph<adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, float>>, typename G::object_type> TempGraph;
	TempGraph f;

#ifdef VERBOSE
	printf("zel input: ");
	for (const auto &p : points) {
		printf("%d,%d ", p.x, p.y);
/*		f.add_vertex(p);*/
	}
	printf("\n");
#endif

	/* construct complete graph */
	for (auto p1 = points.begin(); p1 != points.end(); ++p1) {
		for (auto p2 = std::next(p1, 1); p2 != points.end(); ++p2) {
			typename TempGraph::edge_descriptor e;
			e = f.add_edge(*p1, *p2);
			put(edge_weight, f, e, all_dist[g.vertex(*p1)][g.vertex(*p2)]);
		}
	}
	int n = points.size();
	assert(num_vertices(f) == n && num_edges(f) == n*(n-1)/2);

	/* create all possible triples */
	vector<Triple<TempGraph, G> *> triples;
	int num_combinations = std::max(n*(n-1)*(n-2)/6, 0);
	triples.reserve(num_combinations);
	int comb = 0;

	boost::timer::cpu_timer timer;
	timer.start();

	for (int i = 0; i < n; ++i) {
		for (int j = i+1; j < n; ++j) {
			for (int k = j+1; k < n; ++k) {
				Triple<TempGraph, G> *t = new Triple<TempGraph, G>();
				t->vertices[0] = f.vertex(points[i]);
				t->vertices[1] = f.vertex(points[j]);
				t->vertices[2] = f.vertex(points[k]);

#ifdef VERBOSE
				printf("i[%d]=(%d,%d)[%lu], k[%d]=(%d,%d)[%lu], j[%d]=(%d,%d)[%lu]\n", 
						i, points[i].x, points[i].y, t->vertices[0],
						j, points[j].x, points[j].y, t->vertices[1],
						k, points[k].x, points[k].y, t->vertices[2]);
#endif

				triples.emplace_back(t);
/*				push_heap(triples.begin(), triples.end(), triple_comp<Triple<TempGraph, G> *>());*/
				comb++;
			}
		}
	}
	assert(comb == num_combinations);

	tbb::parallel_for(tbb::blocked_range<int>(0, triples.size(), 128), [&g,&f,&points,&triples,this](const tbb::blocked_range<int> &range) -> void {
		for (int i = range.begin(); i != range.end(); ++i) {
			float min_d = FLT_MAX;
			typename G::vertex_descriptor min_v = G::null_vertex();

			for (const auto &v : g.vertices()) {
				float total_d = 0;
				for (int j = 0; j < 3; ++j) {
					total_d += all_dist[g.vertex(f.object(triples[i]->vertices[j]))][v];
				}
#ifdef VERBOSE
				printf("min: %f current: %f\n", min_d, total_d);
				for (int lol = 0; lol < 3; ++lol) {
					printf("indices: %d ", indices[lol]);
				}
				printf("\n");
#endif
				if (total_d < min_d) {
					min_d = total_d;
					min_v = v;
				}
			}

			assert(min_v != num_vertices(g));
			triples[i]->best_d = min_d;
			triples[i]->best_steiner = min_v;
#ifdef VERBOSE
			printf("best: %f (%d,%d)\n", t->best_d, t->best_steiner.x, t->best_steiner.y);
#endif
		}
	});

	timer.stop();
	boost::timer::cpu_times triple_gen_elapsed = timer.elapsed();
	boost::timer::nanosecond_type triple_gen_time = triple_gen_elapsed.system + triple_gen_elapsed.user;

	/* main algo loop */
	bool done = false;
	set<typename G::object_type> steiners;
	vector<typename TempGraph::vertex_descriptor> mst_pred(num_vertices(f));
	auto mst_pred_map = make_iterator_property_map(mst_pred.begin(), get(vertex_index, f));
	map<typename TempGraph::edge_descriptor, float> contracted_edges;

	timer.start();
	while (!done) {
		const Triple<TempGraph, G> *best_triple = nullptr;
		float max_win = 0;

		prim_minimum_spanning_tree(f, mst_pred_map, weight_map(get(edge_weight, f)));
/*		mst_edge_filter<PredMap> filter(mst_pred_map);*/
/*		filtered_graph<TempGraph, mst_edge_filter<PredMap>> */
		float mst_f_total_before = 0;
		map<typename TempGraph::edge_descriptor, float> mst_before_edges;
#ifdef PRINT_MST_BEFORE 
		printf("MST Before\n");
#endif
		for (const auto &v : f.vertices()) {
			typename TempGraph::vertex_descriptor pre = mst_pred_map[v];
			if (pre != v) {
#ifdef PRINT_MST_BEFORE
			printf("MST edge (%d,%d) -> (%d, %d)\n", f.object(pre).x, f.object(pre).y, f.object(v).x, f.object(v).y);
#endif
				float w = get(edge_weight, f, f.edge(pre, v));
				mst_f_total_before += w;
				assert(mst_before_edges.find(f.edge(pre, v)) == mst_before_edges.end());
				mst_before_edges[f.edge(pre, v)] = w;
			}
		}

		int pruned = 0;

		tbb::spin_mutex g_mutex;
		
/*		for (const auto &t : triples) {*/
		tbb::parallel_for(tbb::blocked_range<int>(0, triples.size(), 1024), [&f,&mst_f_total_before,&triples,&best_triple,&max_win,&g_mutex](const tbb::blocked_range<int> &range) -> void {
/*			printf("f:%X mst_f_total_before:%X triples:%X best_triple:%X max_win:%X g_mutex:%X\n", &f, &mst_f_total_before, &triples, &best_triple, &max_win, &g_mutex);*/
/*			TempGraph *f_local_ptr = new TempGraph(f);*/
/*			TempGraph &f_local = *f_local_ptr;*/
			const TempGraph &f_local = f;
/*			const TempGraph f_local(f);*/
/*			vector<typename TempGraph::vertex_descriptor> *mst_pred_local_ptr = new vector<typename TempGraph::vertex_descriptor>(num_vertices(f_local));*/
/*			vector<typename TempGraph::vertex_descriptor> &mst_pred_local = *mst_pred_local_ptr;*/
			vector<typename TempGraph::vertex_descriptor> mst_pred_local(num_vertices(f_local));
			auto mst_pred_map_local = make_iterator_property_map(mst_pred_local.begin(), get(vertex_index, f_local));

/*			printf("range begin: %d range end: %d\n", range.begin(), range.end());*/

			for (int i = range.begin(); i != range.end(); ++i) {

				const Triple<TempGraph, G> *t = triples[i];

				typename TempGraph::vertex_descriptor u = t->vertices[0];	
				typename TempGraph::vertex_descriptor v = t->vertices[1];	
				typename TempGraph::vertex_descriptor w = t->vertices[2];	

				typename TempGraph::edge_descriptor uv = f_local.edge(u, v);
				typename TempGraph::edge_descriptor vw = f_local.edge(v, w);

/*				float old_uv_length = get(edge_weight, f_local, uv);*/
/*				float old_vw_length = get(edge_weight, f_local, vw);*/

				/* A serious bug in the previous commit caused a large amount of triples to be pruned */
				/* Bug: after pruning a triple, the heap is not updated. therefore, front always remain the same */

#ifdef VERBOSE
				printf("Contracting (%d,%d)[%lu] -> (%d,%d)[%lu] -> (%d,%d)[%lu] old_uv: %f old_vw: %f\n", 
						f_local.object(u).x, f_local.object(u).y, u,
						f_local.object(v).x, f_local.object(v).y, v,
						f_local.object(w).x, f_local.object(w).y, w,
						old_uv_length, old_vw_length);
#endif
				auto filter = edge_weight_filter<TempGraph>(get(edge_weight, f_local), uv, vw);
				prim_minimum_spanning_tree(f_local, mst_pred_map_local, weight_map(filter));

#ifdef PRINT_MST_AFTER 
				printf("MST After\n");
#endif
				typename TempGraph::vertex_iterator vert;
				float mst_f_total_after = 0;
				for (const auto &v : f_local.vertices()) {
					typename TempGraph::vertex_descriptor pre = mst_pred_map_local[v];
					if (pre != v) {
#ifdef PRINT_MST_AFTER 
						printf("MST edge (%d,%d) -> (%d, %d)\n", f_local.object(pre).x, f_local.object(pre).y, f_local.object(v).x, f_local.object(v).y);
#endif
						mst_f_total_after += get(filter, f_local.edge(pre, v));
					}
				}

				float win = mst_f_total_before - mst_f_total_after - t->best_d;

				{
					tbb::spin_mutex::scoped_lock lock(g_mutex);
	/*				g_mutex.lock();*/
					if (win > max_win) {
						best_triple = t;
						max_win = win;
				}
/*				g_mutex.unlock();*/
				}
				/*				float offset = 0;*/
				/*				auto c_uv = contracted_edges.find(uv);*/
				/*				if (c_uv != contracted_edges.end()) {*/
				/*					offset += -c_uv->second;*/
				/*				}*/
				/*				auto c_vw = contracted_edges.find(vw);*/
				/*				if (c_vw != contracted_edges.end()) {*/
				/*					offset += -c_vw->second;*/
				/*				}*/
				/*			map<TempGraph::edge_descriptor, float>::iterator b_uv = mst_before_edges.find(uv);*/
				/*			if (b_uv != mst_before_edges.end()) {*/
				/*				offset = b_uv->second;*/
				/*			}*/
				/*			map<TempGraph::edge_descriptor, float>::iterator b_vw = mst_before_edges.find(vw);*/
				/*			if (b_vw != mst_before_edges.end()) {*/
				/*				offset = b_vw->second;*/
				/*			}*/

#ifdef PRINT_STATS
				printf("mst_before: %f mst_after: %f offset: %f old_uv: %f old_vw: %f\n", mst_f_total_before, mst_f_total_after, offset, old_uv_length, old_vw_length);
#endif

				/*always true*/
				/*				if (mst_f_total_before - mst_f_total_after > old_uv_length+old_vw_length + offset) {*/
				/*					printf("mst_before: %f mst_after: %f old_uv: %f old_vw: %f offset: %f c_uv: %f c_vw: %f\n", mst_f_total_before, mst_f_total_after, old_uv_length, old_vw_length, offset, c_uv->second, c_vw->second);*/
				/*				}*/
				/*				assert(mst_f_total_before - mst_f_total_after <= old_uv_length+old_vw_length + offset); */

				/* false */
				/*			assert(mst_f_total_before-mst_f_total_after >= old_uv_length || mst_f_total_before-mst_f_total_after >= old_vw_length);*/
				/*			assert(mst_f_total_before - mst_f_total_after + offset == old_uv_length+old_vw_length);*/
				/*			assert(mst_f_total_before - mst_f_total_after == old_uv_length+old_vw_length || */
				/*					mst_f_total_before - mst_f_total_after == old_uv_length ||*/
				/*					mst_f_total_before - mst_f_total_after == old_vw_length );*/
#ifdef VERBOSE
				printf("Total before: %f Total after: %f Best d: %f ", mst_f_total_before, mst_f_total_after, t->best_d);
				printf("Win: %f Max win: %f ", win, max_win);
				if (old_uv_length + old_vw_length < t->best_d) {
					pruned++;
					if (win > 0) {
						printf("WEIRD! win=%f when bound is satisfied\n", win);
					}
					/*					assert(win <= 0);*/
					printf("PRUNED");
				}
				printf("\n");

				if (mst_f_total_before - mst_f_total_after > old_uv_length+old_vw_length) {
					printf("BOUND VIOLATION=%f: mst_before: %f mst_after: %f old_uv: %f old_vw: %f\n", mst_f_total_before-mst_f_total_after-(old_uv_length+old_vw_length), mst_f_total_before, mst_f_total_after, old_uv_length, old_vw_length);
				}
				/*				assert(!(mst_f_total_before - mst_f_total_after > old_uv_length + old_vw_length)); */

#endif
				/*			printf("Win: %f best_d: %f\n", win, t->best_d);*/
			}
/*			delete mst_pred_local_ptr;*/
/*			delete f_local_ptr;*/
		});

#ifdef VERBOSE
		printf("Pruned: %d num_triples: %d\n", pruned, triples.size());
#endif

		if (!best_triple) {
			done = true;
		} else {
			typename TempGraph::vertex_descriptor u = best_triple->vertices[0];
			typename TempGraph::vertex_descriptor v = best_triple->vertices[1];
			typename TempGraph::vertex_descriptor w = best_triple->vertices[2];

			typename TempGraph::edge_descriptor uv = f.edge(u, v);
			typename TempGraph::edge_descriptor vw = f.edge(v, w);

			assert(contracted_edges.find(uv) == contracted_edges.end() && contracted_edges.find(vw) == contracted_edges.end());

			contracted_edges[uv] = get(edge_weight, f, uv);
			contracted_edges[vw] = get(edge_weight, f, uv);

			put(edge_weight, f, uv, 0);
			put(edge_weight, f, vw, 0);


/*			printf("Best steiner: (%d,%d)\n", best_triple->best_steiner.x, best_triple->best_steiner.y);*/

			steiners.insert(g.object(best_triple->best_steiner));
		}	
	}
	
	timer.stop();
	boost::timer::cpu_times main_elapsed = timer.elapsed();
	boost::timer::nanosecond_type main_time = main_elapsed.system + main_elapsed.user;

/*	printf("ended\n");*/

#ifdef VERBOSE
	printf("Steiners[count=%d]: ", steiners.size());
	for (const auto &s : steiners) {
		printf("%d,%d ", s.x, s.y);
	}
	printf("\n");
#endif

	if (!previous_steiners.empty()) {
#ifdef VERBOSE
		printf("Previous steiners[count=%d]: ", previous_steiners.size());
		for (const auto &prev_s : previous_steiners) {
			printf("%d,%d ", prev_s.x, prev_s.y);
		}
		printf("\n");
#endif
		for (const auto &s : steiners) {
			float min_d = FLT_MAX;
			const typename G::object_type *min_v = nullptr;
			for (const auto &prev_s : previous_steiners) {
				float d = all_dist[g.vertex(s)][g.vertex(prev_s)];
				if (d < min_d) {
					min_d = d;
					min_v = &prev_s;
				}
			}
#ifdef VERBOSE
			cout << s << " is closest to " << *min_v << " dist=" << min_d << endl;
#endif
		}
		previous_steiners.clear();
	}
/*	for (const auto &s : steiners) {*/
/*		previous_steiners.push_back(s);*/
/*	}*/
/*	printf("done\n");*/

	std::copy(steiners.begin(), steiners.end(), std::back_inserter(previous_steiners));

/*	steiners.insert(points.begin(), points.end());*/

/*	timer.start();*/
/*	list<typename G::edge_descriptor> tree_edges = kmb(g, steiners);*/
/*	timer.stop();*/
/**/
/*	boost::timer::cpu_times kmb_elapsed = timer.elapsed();*/
/*	boost::timer::nanosecond_type kmb_time = kmb_elapsed.system + kmb_elapsed.user;*/
/**/
/*#ifdef VERBOSE*/
/*	printf("Triple gen: %fms Main: %fms Kmb: %fms\n", (double)triple_gen_time/1000000, (double)main_time/1000000, (double)kmb_time/1000000);*/
/*#endif*/

/*	if (!connected) {*/
/*		ofstream file("bug.viz");*/
/*		write_graphviz(file, verify, vertex_label_writer<VerifyGraph>(&verify));*/
/*		exit(0);*/
/*	}*/

	for (const auto &t : triples) {
		delete t;
	}
}

template<typename G>
template<typename Container, typename PreviousSteiner>
void Mst<G>::parallel_zel(const G &g, const Container &points, PreviousSteiner &previous_steiners) const
{
/*	if (do_init) {*/
/*		init(g);*/
/*	}*/

/*	VertexIterator v1, v1_end;*/
/*	VertexIterator v2, v2_end;*/

/*	for (tie(v1, v1_end) = vertices(g); v != v1_end; ++v1) {*/
/*		for (tie(v2, v2_end) = vertices(g); v2 != v2_end; ++v2) {*/
/*			Point p1 = g.object(*v);*/
/*			Point p2 = g.object(*v2);*/
/*			assert(all_dist[*v][*v2] == abs(p1.x - p2.x) + abs(p1.y - p2.y));*/
/*		}*/
/*	}*/

	typedef Graph<adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, float>>, typename G::object_type> TempGraph;
	TempGraph f;

#ifdef VERBOSE
	printf("zel input: ");
	for (const auto &p : points) {
		printf("%d,%d ", p.x, p.y);
/*		f.add_vertex(p);*/
	}
	printf("\n");
#endif

	/* construct complete graph */
	for (auto p1 = points.begin(); p1 != points.end(); ++p1) {
		for (auto p2 = std::next(p1, 1); p2 != points.end(); ++p2) {
			typename TempGraph::edge_descriptor e;
			e = f.add_edge(*p1, *p2);
			put(edge_weight, f, e, all_dist[g.vertex(*p1)][g.vertex(*p2)]);
		}
	}
	int n = points.size();
	assert(num_vertices(f) == n && num_edges(f) == n*(n-1)/2);

	/* create all possible triples */
	vector<Triple<TempGraph, G> *> triples;
	int num_combinations = std::max(n*(n-1)*(n-2)/6, 0);
	triples.reserve(num_combinations);
	int comb = 0;

	boost::timer::cpu_timer timer;
	timer.start();

	for (int i = 0; i < n; ++i) {
		for (int j = i+1; j < n; ++j) {
			for (int k = j+1; k < n; ++k) {
				Triple<TempGraph, G> *t = new Triple<TempGraph, G>();
				t->vertices[0] = f.vertex(points[i]);
				t->vertices[1] = f.vertex(points[j]);
				t->vertices[2] = f.vertex(points[k]);

#ifdef VERBOSE
				printf("i[%d]=(%d,%d)[%lu], k[%d]=(%d,%d)[%lu], j[%d]=(%d,%d)[%lu]\n", 
						i, points[i].x, points[i].y, t->vertices[0],
						j, points[j].x, points[j].y, t->vertices[1],
						k, points[k].x, points[k].y, t->vertices[2]);
#endif

				float min_d = FLT_MAX;
				typename G::vertex_descriptor min_v = num_vertices(g);

				vector<int> indices = {i, j, k};

				for (const auto &v : g.vertices()) {
					float total_d = 0;

					for (const auto &l : indices) {
#ifdef VERBOSE
						printf("[%d](%d,%d) -> (%d,%d) = %f\n", l, points[l].x, points[l].y, g.object(v).x, g.object(v).y, all_dist[g.vertex(points[l])][v]);
#endif
						total_d += all_dist[g.vertex(points[l])][v];
					}
#ifdef VERBOSE
					printf("min: %f current: %f\n", min_d, total_d);
					for (int lol = 0; lol < 3; ++lol) {
						printf("indices: %d ", indices[lol]);
					}
					printf("\n");
#endif
					if (total_d < min_d) {
						min_d = total_d;
						min_v = v;
					}
				}

				assert(min_v != num_vertices(g));
				t->best_d = min_d;
				t->best_steiner = min_v;

#ifdef VERBOSE
				printf("best: %f (%d,%d)\n", t->best_d, t->best_steiner.x, t->best_steiner.y);
#endif

				triples.emplace_back(t);
				push_heap(triples.begin(), triples.end(), triple_comp<Triple<TempGraph, G> *>());
				comb++;
			}
		}
	}
	assert(comb == num_combinations);
	timer.stop();
	boost::timer::cpu_times triple_gen_elapsed = timer.elapsed();
	boost::timer::nanosecond_type triple_gen_time = triple_gen_elapsed.system + triple_gen_elapsed.user;

	/* main algo loop */
	bool done = false;
	set<typename G::object_type> steiners;
	vector<typename TempGraph::vertex_descriptor> mst_pred(num_vertices(f));
	auto mst_pred_map = make_iterator_property_map(mst_pred.begin(), get(vertex_index, f));
	map<typename TempGraph::edge_descriptor, float> contracted_edges;

	timer.start();
	while (!done) {
		const Triple<TempGraph, G> *best_triple = nullptr;
		float max_win = 0;

		prim_minimum_spanning_tree(f, mst_pred_map, weight_map(get(edge_weight, f)));
/*		mst_edge_filter<PredMap> filter(mst_pred_map);*/
/*		filtered_graph<TempGraph, mst_edge_filter<PredMap>> */
		float mst_f_total_before = 0;
		map<typename TempGraph::edge_descriptor, float> mst_before_edges;
#ifdef PRINT_MST_BEFORE 
		printf("MST Before\n");
#endif
		for (const auto &v : f.vertices()) {
			typename TempGraph::vertex_descriptor pre = mst_pred_map[v];
			if (pre != v) {
#ifdef PRINT_MST_BEFORE
			printf("MST edge (%d,%d) -> (%d, %d)\n", f.object(pre).x, f.object(pre).y, f.object(v).x, f.object(v).y);
#endif
				float w = get(edge_weight, f, f.edge(pre, v));
				mst_f_total_before += w;
				assert(mst_before_edges.find(f.edge(pre, v)) == mst_before_edges.end());
				mst_before_edges[f.edge(pre, v)] = w;
			}
		}

		int pruned = 0;

		tbb::spin_mutex g_mutex;
		
/*		for (const auto &t : triples) {*/
		tbb::parallel_for(tbb::blocked_range<int>(0, triples.size(), 1024), [&f,&mst_f_total_before,&triples,&best_triple,&max_win,&g_mutex](const tbb::blocked_range<int> &range) -> void {
/*			printf("f:%X mst_f_total_before:%X triples:%X best_triple:%X max_win:%X g_mutex:%X\n", &f, &mst_f_total_before, &triples, &best_triple, &max_win, &g_mutex);*/
/*			TempGraph *f_local_ptr = new TempGraph(f);*/
/*			TempGraph &f_local = *f_local_ptr;*/
			TempGraph f_local(f);
/*			vector<typename TempGraph::vertex_descriptor> *mst_pred_local_ptr = new vector<typename TempGraph::vertex_descriptor>(num_vertices(f_local));*/
/*			vector<typename TempGraph::vertex_descriptor> &mst_pred_local = *mst_pred_local_ptr;*/
			vector<typename TempGraph::vertex_descriptor> mst_pred_local(num_vertices(f_local));
			auto mst_pred_map_local = make_iterator_property_map(mst_pred_local.begin(), get(vertex_index, f_local));

/*			printf("range begin: %d range end: %d\n", range.begin(), range.end());*/

			for (int i = range.begin(); i != range.end(); ++i) {

				const Triple<TempGraph, G> *t = triples[i];

				typename TempGraph::vertex_descriptor u = t->vertices[0];	
				typename TempGraph::vertex_descriptor v = t->vertices[1];	
				typename TempGraph::vertex_descriptor w = t->vertices[2];	

				typename TempGraph::edge_descriptor uv = f_local.edge(u, v);
				typename TempGraph::edge_descriptor vw = f_local.edge(v, w);

				float old_uv_length = get(edge_weight, f_local, uv);
				float old_vw_length = get(edge_weight, f_local, vw);

				/* A serious bug in the previous commit caused a large amount of triples to be pruned */
				/* Bug: after pruning a triple, the heap is not updated. therefore, front always remain the same */

				put(edge_weight, f_local, uv, 0);
				put(edge_weight, f_local, vw, 0);

#ifdef VERBOSE
				printf("Contracting (%d,%d)[%lu] -> (%d,%d)[%lu] -> (%d,%d)[%lu] old_uv: %f old_vw: %f\n", 
						f_local.object(u).x, f_local.object(u).y, u,
						f_local.object(v).x, f_local.object(v).y, v,
						f_local.object(w).x, f_local.object(w).y, w,
						old_uv_length, old_vw_length);
#endif

				prim_minimum_spanning_tree(f_local, mst_pred_map_local, weight_map(get(edge_weight, f_local)));

#ifdef PRINT_MST_AFTER 
				printf("MST After\n");
#endif
				typename TempGraph::vertex_iterator vert;
				float mst_f_total_after = 0;
				for (const auto &v : f_local.vertices()) {
					typename TempGraph::vertex_descriptor pre = mst_pred_map_local[v];
					if (pre != v) {
#ifdef PRINT_MST_AFTER 
						printf("MST edge (%d,%d) -> (%d, %d)\n", f_local.object(pre).x, f_local.object(pre).y, f_local.object(v).x, f_local.object(v).y);
#endif
						mst_f_total_after += get(edge_weight, f_local, f_local.edge(pre, v));
					}
				}

				put(edge_weight, f_local, uv, old_uv_length);
				put(edge_weight, f_local, vw, old_vw_length);

				float win = mst_f_total_before - mst_f_total_after - t->best_d;

				{
					tbb::spin_mutex::scoped_lock lock(g_mutex);
	/*				g_mutex.lock();*/
					if (win > max_win) {
						best_triple = t;
						max_win = win;
				}
/*				g_mutex.unlock();*/
				}
				/*				float offset = 0;*/
				/*				auto c_uv = contracted_edges.find(uv);*/
				/*				if (c_uv != contracted_edges.end()) {*/
				/*					offset += -c_uv->second;*/
				/*				}*/
				/*				auto c_vw = contracted_edges.find(vw);*/
				/*				if (c_vw != contracted_edges.end()) {*/
				/*					offset += -c_vw->second;*/
				/*				}*/
				/*			map<TempGraph::edge_descriptor, float>::iterator b_uv = mst_before_edges.find(uv);*/
				/*			if (b_uv != mst_before_edges.end()) {*/
				/*				offset = b_uv->second;*/
				/*			}*/
				/*			map<TempGraph::edge_descriptor, float>::iterator b_vw = mst_before_edges.find(vw);*/
				/*			if (b_vw != mst_before_edges.end()) {*/
				/*				offset = b_vw->second;*/
				/*			}*/

#ifdef PRINT_STATS
				printf("mst_before: %f mst_after: %f offset: %f old_uv: %f old_vw: %f\n", mst_f_total_before, mst_f_total_after, offset, old_uv_length, old_vw_length);
#endif

				/*always true*/
				/*				if (mst_f_total_before - mst_f_total_after > old_uv_length+old_vw_length + offset) {*/
				/*					printf("mst_before: %f mst_after: %f old_uv: %f old_vw: %f offset: %f c_uv: %f c_vw: %f\n", mst_f_total_before, mst_f_total_after, old_uv_length, old_vw_length, offset, c_uv->second, c_vw->second);*/
				/*				}*/
				/*				assert(mst_f_total_before - mst_f_total_after <= old_uv_length+old_vw_length + offset); */

				/* false */
				/*			assert(mst_f_total_before-mst_f_total_after >= old_uv_length || mst_f_total_before-mst_f_total_after >= old_vw_length);*/
				/*			assert(mst_f_total_before - mst_f_total_after + offset == old_uv_length+old_vw_length);*/
				/*			assert(mst_f_total_before - mst_f_total_after == old_uv_length+old_vw_length || */
				/*					mst_f_total_before - mst_f_total_after == old_uv_length ||*/
				/*					mst_f_total_before - mst_f_total_after == old_vw_length );*/
#ifdef VERBOSE
				printf("Total before: %f Total after: %f Best d: %f ", mst_f_total_before, mst_f_total_after, t->best_d);
				printf("Win: %f Max win: %f ", win, max_win);
				if (old_uv_length + old_vw_length < t->best_d) {
					pruned++;
					if (win > 0) {
						printf("WEIRD! win=%f when bound is satisfied\n", win);
					}
					/*					assert(win <= 0);*/
					printf("PRUNED");
				}
				printf("\n");

				if (mst_f_total_before - mst_f_total_after > old_uv_length+old_vw_length) {
					printf("BOUND VIOLATION=%f: mst_before: %f mst_after: %f old_uv: %f old_vw: %f\n", mst_f_total_before-mst_f_total_after-(old_uv_length+old_vw_length), mst_f_total_before, mst_f_total_after, old_uv_length, old_vw_length);
				}
				/*				assert(!(mst_f_total_before - mst_f_total_after > old_uv_length + old_vw_length)); */

#endif
				/*			printf("Win: %f best_d: %f\n", win, t->best_d);*/
			}
/*			delete mst_pred_local_ptr;*/
/*			delete f_local_ptr;*/
		});

#ifdef VERBOSE
		printf("Pruned: %d num_triples: %d\n", pruned, triples.size());
#endif

		if (!best_triple) {
			done = true;
		} else {
			typename TempGraph::vertex_descriptor u = best_triple->vertices[0];
			typename TempGraph::vertex_descriptor v = best_triple->vertices[1];
			typename TempGraph::vertex_descriptor w = best_triple->vertices[2];

			typename TempGraph::edge_descriptor uv = f.edge(u, v);
			typename TempGraph::edge_descriptor vw = f.edge(v, w);

			assert(contracted_edges.find(uv) == contracted_edges.end() && contracted_edges.find(vw) == contracted_edges.end());

			contracted_edges[uv] = get(edge_weight, f, uv);
			contracted_edges[vw] = get(edge_weight, f, uv);

			put(edge_weight, f, uv, 0);
			put(edge_weight, f, vw, 0);


/*			printf("Best steiner: (%d,%d)\n", best_triple->best_steiner.x, best_triple->best_steiner.y);*/

			steiners.insert(g.object(best_triple->best_steiner));
		}	
	}
	
	timer.stop();
	boost::timer::cpu_times main_elapsed = timer.elapsed();
	boost::timer::nanosecond_type main_time = main_elapsed.system + main_elapsed.user;

/*	printf("ended\n");*/

#ifdef VERBOSE
	printf("Steiners[count=%d]: ", steiners.size());
	for (const auto &s : steiners) {
		printf("%d,%d ", s.x, s.y);
	}
	printf("\n");
#endif

	if (!previous_steiners.empty()) {
#ifdef VERBOSE
		printf("Previous steiners[count=%d]: ", previous_steiners.size());
		for (const auto &prev_s : previous_steiners) {
			printf("%d,%d ", prev_s.x, prev_s.y);
		}
		printf("\n");
#endif
		for (const auto &s : steiners) {
			float min_d = FLT_MAX;
			const typename G::object_type *min_v = nullptr;
			for (const auto &prev_s : previous_steiners) {
				float d = all_dist[g.vertex(s)][g.vertex(prev_s)];
				if (d < min_d) {
					min_d = d;
					min_v = &prev_s;
				}
			}
#ifdef VERBOSE
			cout << s << " is closest to " << *min_v << " dist=" << min_d << endl;
#endif
		}
		previous_steiners.clear();
	}
/*	for (const auto &s : steiners) {*/
/*		previous_steiners.push_back(s);*/
/*	}*/
/*	printf("done\n");*/

	std::copy(steiners.begin(), steiners.end(), std::back_inserter(previous_steiners));

/*	steiners.insert(points.begin(), points.end());*/

/*	timer.start();*/
/*	list<typename G::edge_descriptor> tree_edges = kmb(g, steiners);*/
/*	timer.stop();*/
/**/
/*	boost::timer::cpu_times kmb_elapsed = timer.elapsed();*/
/*	boost::timer::nanosecond_type kmb_time = kmb_elapsed.system + kmb_elapsed.user;*/
/**/
/*#ifdef VERBOSE*/
/*	printf("Triple gen: %fms Main: %fms Kmb: %fms\n", (double)triple_gen_time/1000000, (double)main_time/1000000, (double)kmb_time/1000000);*/
/*#endif*/

/*	if (!connected) {*/
/*		ofstream file("bug.viz");*/
/*		write_graphviz(file, verify, vertex_label_writer<VerifyGraph>(&verify));*/
/*		exit(0);*/
/*	}*/

	for (const auto &t : triples) {
		delete t;
	}
}

template<typename G>
template<typename Container, typename PreviousSteiner>
list<typename G::edge_descriptor> Mst<G>::zel(const G &g, const Container &points, PreviousSteiner &previous_steiners) const
{
/*	if (do_init) {*/
/*		init(g);*/
/*	}*/

/*	VertexIterator v1, v1_end;*/
/*	VertexIterator v2, v2_end;*/

/*	for (tie(v1, v1_end) = vertices(g); v != v1_end; ++v1) {*/
/*		for (tie(v2, v2_end) = vertices(g); v2 != v2_end; ++v2) {*/
/*			Point p1 = g.object(*v);*/
/*			Point p2 = g.object(*v2);*/
/*			assert(all_dist[*v][*v2] == abs(p1.x - p2.x) + abs(p1.y - p2.y));*/
/*		}*/
/*	}*/

	typedef Graph<adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, float>>, typename G::object_type> TempGraph;
	TempGraph f;

#ifdef VERBOSE
	printf("zel input: ");
	for (const auto &p : points) {
		printf("%d,%d ", p.x, p.y);
/*		f.add_vertex(p);*/
	}
	printf("\n");
#endif

	/* construct complete graph */
	for (auto p1 = points.begin(); p1 != points.end(); ++p1) {
		for (auto p2 = std::next(p1, 1); p2 != points.end(); ++p2) {
			typename TempGraph::edge_descriptor e;
			e = f.add_edge(*p1, *p2);
			put(edge_weight, f, e, all_dist[g.vertex(*p1)][g.vertex(*p2)]);
		}
	}
	int n = points.size();
	assert(num_vertices(f) == n && num_edges(f) == n*(n-1)/2);

	/* create all possible triples */
	vector<Triple<TempGraph, G> *> triples;
	int num_combinations = std::max(n*(n-1)*(n-2)/6, 0);
	triples.reserve(num_combinations);
	int comb = 0;

	boost::timer::cpu_timer timer;
	timer.start();

	for (int i = 0; i < n; ++i) {
		for (int j = i+1; j < n; ++j) {
			for (int k = j+1; k < n; ++k) {
				Triple<TempGraph, G> *t = new Triple<TempGraph, G>();
				t->vertices[0] = f.vertex(points[i]);
				t->vertices[1] = f.vertex(points[j]);
				t->vertices[2] = f.vertex(points[k]);

#ifdef VERBOSE
				printf("i[%d]=(%d,%d)[%lu], k[%d]=(%d,%d)[%lu], j[%d]=(%d,%d)[%lu]\n", 
						i, points[i].x, points[i].y, t->vertices[0],
						j, points[j].x, points[j].y, t->vertices[1],
						k, points[k].x, points[k].y, t->vertices[2]);
#endif

				float min_d = FLT_MAX;
				typename G::vertex_descriptor min_v = num_vertices(g);

				vector<int> indices = {i, j, k};

				for (const auto &v : g.vertices()) {
					float total_d = 0;

					for (const auto &l : indices) {
#ifdef VERBOSE
						printf("[%d](%d,%d) -> (%d,%d) = %f\n", l, points[l].x, points[l].y, g.object(v).x, g.object(v).y, all_dist[g.vertex(points[l])][v]);
#endif
						total_d += all_dist[g.vertex(points[l])][v];
					}
#ifdef VERBOSE
					printf("min: %f current: %f\n", min_d, total_d);
					for (int lol = 0; lol < 3; ++lol) {
						printf("indices: %d ", indices[lol]);
					}
					printf("\n");
#endif
					if (total_d < min_d) {
						min_d = total_d;
						min_v = v;
					}
				}

				assert(min_v != num_vertices(g));
				t->best_d = min_d;
				t->best_steiner = min_v;

#ifdef VERBOSE
				printf("best: %f (%d,%d)\n", t->best_d, t->best_steiner.x, t->best_steiner.y);
#endif

				triples.emplace_back(t);
				push_heap(triples.begin(), triples.end(), triple_comp<Triple<TempGraph, G> *>());
				comb++;
			}
		}
	}
	assert(comb == num_combinations);
	timer.stop();
	boost::timer::cpu_times triple_gen_elapsed = timer.elapsed();
	boost::timer::nanosecond_type triple_gen_time = triple_gen_elapsed.system + triple_gen_elapsed.user;

	/* main algo loop */
	bool done = false;
	set<typename G::object_type> steiners;
	vector<typename TempGraph::vertex_descriptor> mst_pred(num_vertices(f));
	auto mst_pred_map = make_iterator_property_map(mst_pred.begin(), get(vertex_index, f));
	map<typename TempGraph::edge_descriptor, float> contracted_edges;

	timer.start();
	while (!done) {
		const Triple<TempGraph, G> *best_triple = nullptr;
		float max_win = 0;

		prim_minimum_spanning_tree(f, mst_pred_map, weight_map(get(edge_weight, f)));
/*		mst_edge_filter<PredMap> filter(mst_pred_map);*/
/*		filtered_graph<TempGraph, mst_edge_filter<PredMap>> */
		float mst_f_total_before = 0;
		map<typename TempGraph::edge_descriptor, float> mst_before_edges;
#ifdef PRINT_MST_BEFORE 
		printf("MST Before\n");
#endif
		for (const auto &v : f.vertices()) {
			typename TempGraph::vertex_descriptor pre = mst_pred_map[v];
			if (pre != v) {
#ifdef PRINT_MST_BEFORE
			printf("MST edge (%d,%d) -> (%d, %d)\n", f.object(pre).x, f.object(pre).y, f.object(v).x, f.object(v).y);
#endif
				float w = get(edge_weight, f, f.edge(pre, v));
				mst_f_total_before += w;
				assert(mst_before_edges.find(f.edge(pre, v)) == mst_before_edges.end());
				mst_before_edges[f.edge(pre, v)] = w;
			}
		}

		int pruned = 0;
		
		make_heap(triples.begin(), triples.end(), triple_comp<Triple<TempGraph, G> *>());

/*		for (const auto &t : triples) {*/
		for (int i = 0; i < triples.size(); ++i) {
			const Triple<TempGraph, G> *t = triples.front();

			typename TempGraph::vertex_descriptor u = t->vertices[0];	
			typename TempGraph::vertex_descriptor v = t->vertices[1];	
			typename TempGraph::vertex_descriptor w = t->vertices[2];	

			typename TempGraph::edge_descriptor uv = f.edge(u, v);
			typename TempGraph::edge_descriptor vw = f.edge(v, w);

			float old_uv_length = get(edge_weight, f, uv);
			float old_vw_length = get(edge_weight, f, vw);
			
			/* A serious bug in the previous commit caused a large amount of triples to be pruned */
			/* Bug: after pruning a triple, the heap is not updated. therefore, front always remain the same */
			/*if (old_uv_length + old_vw_length < t->best_d) {
				printf("Total before: %f Total after: %f Best d: %f ", mst_f_total_before, mst_f_total_after, t->best_d);
				pruned++;
			} else */{

				put(edge_weight, f, uv, 0);
				put(edge_weight, f, vw, 0);

#ifdef VERBOSE
				printf("Contracting (%d,%d)[%lu] -> (%d,%d)[%lu] -> (%d,%d)[%lu] old_uv: %f old_vw: %f\n", 
						f.object(u).x, f.object(u).y, u,
						f.object(v).x, f.object(v).y, v,
						f.object(w).x, f.object(w).y, w,
						old_uv_length, old_vw_length);
#endif

				prim_minimum_spanning_tree(f, mst_pred_map, weight_map(get(edge_weight, f)));

#ifdef PRINT_MST_AFTER 
				printf("MST After\n");
#endif
				typename TempGraph::vertex_iterator vert;
				float mst_f_total_after = 0;
				for (const auto &v : f.vertices()) {
					typename TempGraph::vertex_descriptor pre = mst_pred_map[v];
					if (pre != v) {
#ifdef PRINT_MST_AFTER 
						printf("MST edge (%d,%d) -> (%d, %d)\n", f.object(pre).x, f.object(pre).y, f.object(v).x, f.object(v).y);
#endif
						mst_f_total_after += get(edge_weight, f, f.edge(pre, v));
					}
				}

				put(edge_weight, f, uv, old_uv_length);
				put(edge_weight, f, vw, old_vw_length);

				float win = mst_f_total_before - mst_f_total_after - t->best_d;
				if (win > max_win) {
					best_triple = t;
					max_win = win;
				}
/*				float offset = 0;*/
/*				auto c_uv = contracted_edges.find(uv);*/
/*				if (c_uv != contracted_edges.end()) {*/
/*					offset += -c_uv->second;*/
/*				}*/
/*				auto c_vw = contracted_edges.find(vw);*/
/*				if (c_vw != contracted_edges.end()) {*/
/*					offset += -c_vw->second;*/
/*				}*/
	/*			map<TempGraph::edge_descriptor, float>::iterator b_uv = mst_before_edges.find(uv);*/
	/*			if (b_uv != mst_before_edges.end()) {*/
	/*				offset = b_uv->second;*/
	/*			}*/
	/*			map<TempGraph::edge_descriptor, float>::iterator b_vw = mst_before_edges.find(vw);*/
	/*			if (b_vw != mst_before_edges.end()) {*/
	/*				offset = b_vw->second;*/
	/*			}*/
					
#ifdef PRINT_STATS
				printf("mst_before: %f mst_after: %f offset: %f old_uv: %f old_vw: %f\n", mst_f_total_before, mst_f_total_after, offset, old_uv_length, old_vw_length);
#endif

				/*always true*/
/*				if (mst_f_total_before - mst_f_total_after > old_uv_length+old_vw_length + offset) {*/
/*					printf("mst_before: %f mst_after: %f old_uv: %f old_vw: %f offset: %f c_uv: %f c_vw: %f\n", mst_f_total_before, mst_f_total_after, old_uv_length, old_vw_length, offset, c_uv->second, c_vw->second);*/
/*				}*/
/*				assert(mst_f_total_before - mst_f_total_after <= old_uv_length+old_vw_length + offset); */

				/* false */
	/*			assert(mst_f_total_before-mst_f_total_after >= old_uv_length || mst_f_total_before-mst_f_total_after >= old_vw_length);*/
	/*			assert(mst_f_total_before - mst_f_total_after + offset == old_uv_length+old_vw_length);*/
	/*			assert(mst_f_total_before - mst_f_total_after == old_uv_length+old_vw_length || */
	/*					mst_f_total_before - mst_f_total_after == old_uv_length ||*/
	/*					mst_f_total_before - mst_f_total_after == old_vw_length );*/
#ifdef VERBOSE
				printf("Total before: %f Total after: %f Best d: %f ", mst_f_total_before, mst_f_total_after, t->best_d);
				printf("Win: %f Max win: %f ", win, max_win);
				if (old_uv_length + old_vw_length < t->best_d) {
					pruned++;
					if (win > 0) {
						printf("WEIRD! win=%f when bound is satisfied\n", win);
					}
/*					assert(win <= 0);*/
					printf("PRUNED");
				}
				printf("\n");

				if (mst_f_total_before - mst_f_total_after > old_uv_length+old_vw_length) {
					printf("BOUND VIOLATION=%f: mst_before: %f mst_after: %f old_uv: %f old_vw: %f\n", mst_f_total_before-mst_f_total_after-(old_uv_length+old_vw_length), mst_f_total_before, mst_f_total_after, old_uv_length, old_vw_length);
				}
/*				assert(!(mst_f_total_before - mst_f_total_after > old_uv_length + old_vw_length)); */

#endif
	/*			printf("Win: %f best_d: %f\n", win, t->best_d);*/
			}

			pop_heap(triples.begin(), triples.end()-i, triple_comp<Triple<TempGraph, G> *>());
		}

#ifdef VERBOSE
		printf("Pruned: %d num_triples: %d\n", pruned, triples.size());
#endif

		if (!best_triple) {
			done = true;
		} else {
			typename TempGraph::vertex_descriptor u = best_triple->vertices[0];
			typename TempGraph::vertex_descriptor v = best_triple->vertices[1];
			typename TempGraph::vertex_descriptor w = best_triple->vertices[2];

			typename TempGraph::edge_descriptor uv = f.edge(u, v);
			typename TempGraph::edge_descriptor vw = f.edge(v, w);

			assert(contracted_edges.find(uv) == contracted_edges.end() && contracted_edges.find(vw) == contracted_edges.end());

			contracted_edges[uv] = get(edge_weight, f, uv);
			contracted_edges[vw] = get(edge_weight, f, uv);

			put(edge_weight, f, uv, 0);
			put(edge_weight, f, vw, 0);


/*			printf("Best steiner: (%d,%d)\n", best_triple->best_steiner.x, best_triple->best_steiner.y);*/

			steiners.insert(g.object(best_triple->best_steiner));
		}	
	}
	
	timer.stop();
	boost::timer::cpu_times main_elapsed = timer.elapsed();
	boost::timer::nanosecond_type main_time = main_elapsed.system + main_elapsed.user;


#ifdef VERBOSE
	printf("Steiners[count=%d]: ", steiners.size());
	for (const auto &s : steiners) {
		printf("%d,%d ", s.x, s.y);
	}
	printf("\n");
#endif

	if (!previous_steiners.empty()) {
#ifdef VERBOSE
		printf("Previous steiners[count=%d]: ", previous_steiners.size());
		for (const auto &prev_s : previous_steiners) {
			printf("%d,%d ", prev_s.x, prev_s.y);
		}
		printf("\n");
#endif
		for (const auto &s : steiners) {
			float min_d = FLT_MAX;
			const typename G::object_type *min_v = nullptr;
			for (const auto &prev_s : previous_steiners) {
				float d = all_dist[g.vertex(s)][g.vertex(prev_s)];
				if (d < min_d) {
					min_d = d;
					min_v = &prev_s;
				}
			}
#ifdef VERBOSE
			cout << s << " is closest to " << *min_v << " dist=" << min_d << endl;
#endif
		}
		previous_steiners.clear();
	}
	for (const auto &s : steiners) {
		previous_steiners.push_back(s);
	}
/*	printf("done\n");*/

	steiners.insert(points.begin(), points.end());

	timer.start();
	list<typename G::edge_descriptor> tree_edges = kmb(g, steiners);
	timer.stop();

	boost::timer::cpu_times kmb_elapsed = timer.elapsed();
	boost::timer::nanosecond_type kmb_time = kmb_elapsed.system + kmb_elapsed.user;

#ifdef VERBOSE
	printf("Triple gen: %fms Main: %fms Kmb: %fms\n", (double)triple_gen_time/1000000, (double)main_time/1000000, (double)kmb_time/1000000);
#endif

/*	if (!connected) {*/
/*		ofstream file("bug.viz");*/
/*		write_graphviz(file, verify, vertex_label_writer<VerifyGraph>(&verify));*/
/*		exit(0);*/
/*	}*/

	for (const auto &t : triples) {
		delete t;
	}

	return tree_edges;
}

/*template class Mst<GridGraph>;*/
/*template void Mst<GridGraph>::zel(const GridGraph &g, const vector<Point> &points, bool do_init);*/
/*template void Mst<Graph<adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, float>>, int>>::kmb(const Graph<adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, float>>, int> &g, const vector<int> &points, bool do_init=true);*/
#endif
