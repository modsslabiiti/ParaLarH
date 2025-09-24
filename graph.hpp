#ifndef GRAPH_H
#define GRAPH_H

using namespace std;
using namespace boost;

template<class Graph>
bool is_connected(const Graph &g)
{
	vector<int> comp(num_vertices(g));
	auto comp_map = make_iterator_property_map(comp.begin(), get(vertex_index, g));
	return connected_components(g, comp_map) == 1;
/*	bool connected = true;*/
/*	for (const auto &c : comp) {*/
/*		if (c != 0) {*/
/*			connected = false;*/
/**/
/*		}*/
/*	}*/
/*	return connected;*/
}

template<class Graph>
bool is_tree(const Graph &g)
{
	return is_connected(g) && num_vertices(g)-1 == num_edges(g);
}

/*template <class OutEdgeListS = vecS, // a Sequence or an AssociativeContainer*/
/*    class VertexListS = vecS, // a Sequence or a RandomAccessContainer*/
/*    class DirectedS = directedS,*/
/*    class Object = no_property,*/
/*    class EdgeProperty = no_property,*/
/*    class GraphProperty = no_property,*/
/*    class EdgeListS = listS>*/
template<class Base, class Object>
class Graph : public Base {
	public:
/*		typedef adjacency_list<OutEdgeListS, VertexListS, DirectedS, Object, EdgeProperty, GraphProperty, EdgeListS> Base;*/
		typedef typename Base::vertex_descriptor vertex_descriptor;
		typedef typename Base::vertex_iterator vertex_iterator;
		typedef typename Base::edge_descriptor edge_descriptor;
		typedef typename Base::edge_iterator edge_iterator;
		typedef Object object_type;

		template<typename T>
			struct Iterable {
				T _begin;
				T _end;

				Iterable(const T &_begin, const T &_end) :
					_begin(_begin), _end(_end)
				{
				}

				const T &begin() const {
					return _begin;
				}
				const T &end() const {
					return _end;
				}
			};

		vertex_descriptor add_vertex(const Object &o) {
			typename map<Object, vertex_descriptor>::const_iterator iter = o_to_v.find(o);

			vertex_descriptor v = graph_traits<Base>::null_vertex();
			if (iter != o_to_v.end()) {
				v = iter->second;
			} else {
				v = boost::add_vertex(*this);
				o_to_v[o] = v;
				v_to_o[v] = o;
			}
			return v;
		}

		vertex_descriptor vertex(const Object &o) const {
			typename map<Object, vertex_descriptor>::const_iterator iter = o_to_v.find(o);

			vertex_descriptor v = graph_traits<Base>::null_vertex();
			if (iter != o_to_v.end()) {
				v = iter->second;
			}
			return v;
		}
		
		const Object &object(const vertex_descriptor &v) const {
			return v_to_o.at(v);
		}

		edge_descriptor add_edge(const Object &u_o, const Object &v_o) {
			vertex_descriptor u = add_vertex(u_o);
			vertex_descriptor v = add_vertex(v_o);

			edge_descriptor e;
			bool exists;
			std::tie(e, exists) = boost::edge(u, v, *this);
			if (!exists) {
				bool added;
				std::tie(e, added) = boost::add_edge(u, v, *this);
				if (!added) {
					throw;
				}
			} else {
/*				throw;*/
			}

			return e;
		}

		edge_descriptor add_edge(const vertex_descriptor &u, const vertex_descriptor &v) {
			edge_descriptor e;
			bool exists;
			std::tie(e, exists) = boost::edge(u, v, *this);
			if (!exists) {
				bool added;
				std::tie(e, added) = boost::add_edge(u, v, *this);
				if (!added) {
					throw;
				}
			} else {
/*				throw;*/
			}

			return e;
		}

		edge_descriptor edge(const Object &u_o, const Object &v_o) const {
			vertex_descriptor u = o_to_v.at(u_o);
			vertex_descriptor v = o_to_v.at(v_o);

			edge_descriptor e;
			bool exists;
			std::tie(e, exists) = boost::edge(u, v, *this);
			if (!exists) {
				throw;
			}
			return e;
		}

		edge_descriptor edge(const vertex_descriptor &u, const vertex_descriptor &v) const {
			edge_descriptor e;
			bool found;
			std::tie(e, found) = boost::edge(u, v, *this);
			if (!found) {
				throw;
			}
			return e;
		}

		Iterable<vertex_iterator> vertices() const {
			vertex_iterator begin, end;
			std::tie(begin, end) = boost::vertices(*this);

			return Iterable<vertex_iterator>(begin, end);
		}

		Iterable<edge_iterator> edges() const {
			edge_iterator begin, end;
			std::tie(begin, end) = boost::edges(*this);

			return Iterable<edge_iterator>(begin, end);
		}

	private:
		map<Object, vertex_descriptor> o_to_v;
		map<vertex_descriptor, Object> v_to_o;
};

#endif
