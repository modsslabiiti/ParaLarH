#include "pch.hpp"
#include "router.hpp"
#include "netlist.hpp"
#include "mst.hpp"
#include "timing.hpp"
#include <cmath>

GridGraph create_grid_graph(int nx, int ny)
{
	GridGraph g;

	g[graph_bundle].nx = nx;
	g[graph_bundle].ny = ny;

	format fmt("%d,%d");
	for (int y = 0; y < ny; ++y) {
		for (int x = 0; x < nx; ++x) {
			fmt % x % y;
			GridGraph::vertex_descriptor v = g.add_vertex(Point(x, y));
			put(vertex_name, g, v, fmt.str());
/*			printf("vertex: %s\n", get(vertex_name, g, v).c_str());*/
		}
	}
	for (int y = 0; y < ny; ++y) {
		for (int x = 0; x < nx; ++x) {
			GridGraph::edge_descriptor e;
			if (x < nx - 1) {
				e = g.add_edge(Point(x, y), Point(x+1, y));
				put(edge_weight, g, e, 1);
			}
			if (y < ny - 1) {
				e = g.add_edge(Point(x, y), Point(x, y+1));
				put(edge_weight, g, e, 1);
			}
		}
	}
	assert(num_vertices(g) == nx*ny && num_edges(g) == 2*(nx-1)*(ny-1)+(nx-1)+(ny-1));
/*	vector_property_map<int, property_map<GridGraph, vertex_index_t>::type> pmap = make_vector_property_map<int>(get(vertex_index, g));*/
/*	printf("num veritices: %d\n", num_vertices(g));*/
/*	printf("vertex:%d\n", vertex(1000,g));*/
/*	pmap[vertex(1000,g )] = 100;*/
/*	printf("size: %d\n")*/
	return g;
}

string current_net_details;

extern int MAX_ITER;
extern int MAX_UTIL;


// dhaarna - start 
//defining data types 
typedef boost::tuples::tuple<double, double, double, double> Edge; // starting and end point of edge

struct ihash
    : std::unary_function<Edge, std::size_t> 			// to use Edge as a key , a hash function is defined	
{
    std::size_t operator()(Edge const& e) const
    {
        std::size_t seed = 0;
        boost::hash_combine( seed, e.get<0>() );
        boost::hash_combine( seed, e.get<1>() );
        boost::hash_combine( seed, e.get<2>() );
        boost::hash_combine( seed, e.get<3>() );
        return seed;
    }
};

struct iequal_to
    : std::binary_function<Edge, Edge, bool>	   // equality operator for Edge data type
{    bool operator()(Edge const& x, Edge const& y) const
    {
        return ( x.get<0>()==y.get<0>() &&
                 x.get<1>()==y.get<1>() &&
                 x.get<2>()==y.get<2>() &&
                 x.get<3>()==y.get<3>());
    }
};

typedef boost::unordered_map< Edge,set<int>, ihash, iequal_to > EdgeMap;	//to store nets using a particular edge

typedef boost::tuples::tuple<double, double, double, double,double> Edges; // net number , starting and end points of edge	

struct ihashs
    : std::unary_function<Edges, std::size_t> 			// to use Edges as a key , a hash function is defined	
{
    std::size_t operator()(Edges const& e) const
    {
        std::size_t seed = 0;
        boost::hash_combine( seed, e.get<0>() );
        boost::hash_combine( seed, e.get<1>() );
        boost::hash_combine( seed, e.get<2>() );
        boost::hash_combine( seed, e.get<3>() );
	boost::hash_combine( seed, e.get<4>() );
        return seed;
    }
};

struct iequal_tos					// equality operator for Edges data type
    : std::binary_function<Edges, Edges, bool>	   //dhaarna		
{
    bool operator()(Edges const& x, Edges const& y) const
    {
        return ( x.get<0>()==y.get<0>() &&
                 x.get<1>()==y.get<1>() &&
                 x.get<2>()==y.get<2>() &&
                 x.get<3>()==y.get<3>()&&
                 x.get<4>()==y.get<4>());
    }
};

typedef boost::unordered_map< Edges, int, ihashs, iequal_tos > NetMap;		//dhaarna // to check if a given edge is present in a map or not
typedef boost::unordered_map < pair<int,int> , vector<pair<int,int> > > PairMap; // to store neighbours of a point in graph
typedef boost::unordered_map < pair<int,int> , int > VisMap;  // utility data structure for bfs to keep track of visited nodes
typedef boost::unordered_map < pair<int,int> , pair<int,int> > ParMap; // to store parents in path found by bfs

//dhaarna -end 
void Router::parallel_route_tbb_new_netlist(Netlist &netlist, int num_threads)
{int minw=INT_MAX; //dhaarna  // minimum channel width (to be calculated)
PairMap vic; //dhaarna
	Mst<GridGraph> mst;
	int max_iter = MAX_ITER;
	boost::timer::nanosecond_type total_routing_time = 0;

	for (const auto &e : grid.edges()) {
		GridGraphEdgeObject &obj = get(edge_object, grid, e);
		obj.mult = 0;
		obj.max_util = MAX_UTIL;
		obj.delay = 1;
		obj.slp=0;
		obj.m=0;
		obj.T1=0;
		obj.T2=0;

vic[make_pair(grid.object(source(e, grid)).x,grid.object(source(e, grid)).y)].push_back(make_pair(grid.object(target(e, grid)).x,grid.object(target(e, grid)).y)); //dhaarna // assigning neighbours
vic[make_pair(grid.object(target(e, grid)).x,grid.object(target(e, grid)).y)].push_back(make_pair(grid.object(source(e, grid)).x,grid.object(source(e, grid)).y)); //dhaarna

	}

	vector<Net> nets;
	nets.reserve(netlist.nets.size());
	int temp = 0;
	for (const auto &net : netlist.nets) {
		Net temp_net;	
		temp_net.index = temp++;
		temp_net.name = net.second->name;
		temp_net.points.emplace_back(net.second->source->port->block->position);

		for (const auto &sink : net.second->sinks) {
			temp_net.points.emplace_back(sink->port->block->position);
		}
/*		printf("%d\n", temp_net.points.size());*/
		nets.emplace_back(std::move(temp_net));
	}

/*	return;*/

	vector<list<Point>> previous_steiners(nets.size());
	int best_iter = -1;
	int best_util = INT_MAX;
	int best_wirelength = -1;
	float best_crit = -1;
	float best_average = -1;
	float ult_crit = FLT_MAX;
	int ult_crit_iter = -1;
	
	//int flag_2=1;
	//int m=0;
	//int t=0;
	//int t_for_min=0;
	//float slp_residual_array[20]={0.0};
	//int max_util_array[20]={0};
	//float slp_residual_min=99999;			// assign maximun float number
	//float change_value=0.0;
	//float change_value_minimum=1.0;
	//int iter=0;
	float rho=0.0001;
	for (int iter = 0; iter < max_iter; ++iter) {

int maxew=0; //dhaarna // maximum edge weight needed at the end of iteration ( which is used to update minw)
EdgeMap enets;			//dhaarna
NetMap netse; //dhaarna
	//while (flag_2) {
		/* best step size so far, 0.01/(iter+1) */
		//float step_size = (float)0.01/(iter+1);
		//iter=iter+1;

		for (const auto &e : grid.edges()) {
			GridGraphEdgeObject &obj = get(edge_object, grid, e);
			obj.util = 0;
			//obj.util_previous = 0;
			assert(obj.delay == 1);
			put(edge_weight, grid, e, obj.delay + obj.mult);
		}
		mst.init(grid);

		vector<list<GridGraph::edge_descriptor>> route_edges(nets.size());

		boost::timer::cpu_timer timer;
		clock_t iter_begin = clock();
		timer.start();

		int total_wirelength = 0;
		const bool nested_zel = true;
		if (iter == 0) {
			if (!nested_zel) {
				int net_num = 0;

				for (const auto &net : nets) {
					mst.parallel_zel_new_reduce(grid, net.points, previous_steiners[net_num]);
					net_num++;
				}
			} else {
				/* previous grain size: 16 */
				tbb::parallel_for(tbb::blocked_range<int>(0, nets.size()), [&mst, &nets, &previous_steiners, this](const tbb::blocked_range<int> &range) -> void
					{
						for (int i = range.begin(); i != range.end(); ++i) {
							if (nets[i].points.size() <= 70) {
								mst.parallel_zel_new_reduce(grid, nets[i].points, previous_steiners[i]);
							}
						}
					});
			}
		}

		tbb::spin_mutex mutex;
		/* previous grain size: 16 */
		tbb::parallel_for(tbb::blocked_range<int>(0, nets.size()), [&mutex, &route_edges, &total_wirelength, &mst, &nets, &previous_steiners, this](const tbb::blocked_range<int> &range) -> void
			{
				for (int i = range.begin(); i != range.end(); ++i) {
/*					if (nets[i].points.size() <= 30) {*/
/*						mst.parallel_zel_new_reduce(grid, nets[i].points, previous_steiners[i]);*/
/*					}*/
					set<Point> steiners;
					steiners.insert(previous_steiners[i].begin(), previous_steiners[i].end());
					steiners.insert(nets[i].points.begin(), nets[i].points.end());
					route_edges[i] = std::move(mst.kmb(grid, steiners));

/*					set<GridGraph::edge_descriptor> route_edges_check;*/
/*					for (const auto &e : route_edges[i]) {*/
/*						route_edges_check.insert(e);*/
/*					}*/
/*					if (route_edges_check.size() != route_edges[i].size()) {*/
/*						printf("check: %d real: %d\n", route_edges_check.size(), route_edges[i].size());*/
/*					}*/
/*					assert(route_edges_check.size() == route_edges[i].size());*/

					{
						tbb::spin_mutex::scoped_lock lock(mutex);
						total_wirelength += route_edges[i].size();
					}
					for (const auto &e : route_edges[i]) {
						GridGraphEdgeObject &obj = get(edge_object, grid, e);
						{
							tbb::spin_mutex::scoped_lock lock(mutex);
							obj.util += 1;
						}
					}
				}
			}
		);

/***************** calculating f=(c+lambda)*x -w*lambda    *************************
		float f_value=0.0;
		for (const auto &e : grid.edges()) {
			GridGraphEdgeObject &obj = get(edge_object, grid, e);
			f_value=f_value+obj.util*(obj.delay + obj.mult)-MAX_UTIL*obj.mult;
		}
		printf("At iteration number %d, the objective function value is %f\n",iter,f_value);
/***************************************************************************/

/******************step size calculation **********************************/
		float norm_T1T2=0.0;

		for (const auto &e : grid.edges()) {
			GridGraphEdgeObject &obj = get(edge_object, grid, e);
			if (obj.util>obj.max_util){
				obj.m=1;
			}
			else {
				obj.m=0;
			}	
			obj.slp=std::max(0.0f, (float) (obj.util-obj.max_util));
			obj.T1=obj.delay+obj.m*(obj.mult+rho*(obj.slp));
			obj.T2=-obj.slp;
			norm_T1T2=norm_T1T2+sqrt(obj.T1*obj.T1+obj.T2*obj.T2);	
		}
	/*	for (const auto &e : grid.edges()) {
			GridGraphEdgeObject &obj = get(edge_object, grid, e);
			obj.slp=std::max(0.0f, (float) (obj.util-obj.max_util));
			obj.T1=obj.delay+obj.m*(obj.mult+rho*(obj.slp));
		}*/
	/*	for (const auto &e : grid.edges()) {
			GridGraphEdgeObject &obj = get(edge_object, grid, e);
			obj.slp=std::max(0.0f, (float) (obj.util-obj.max_util));
			obj.T2=-obj.slp;
		}*/
	/*	float norm_T1T2=0.0;
		for (const auto &e : grid.edges()) {
			GridGraphEdgeObject &obj = get(edge_object, grid, e);
			norm_T1T2=norm_T1T2+sqrt(obj.T1*obj.T1+obj.T2*obj.T2);
		//printf("for iter=%d m=%d slp=%f T1=%f T2=%f norm of T1;T2=%f\n",iter,obj.m,obj.slp,obj.T1,obj.T2,norm_T1T2);
		}*/
		float step_size = ((float)1.01/(iter+1))/(norm_T1T2);
//		printf("for iter=%d, the step size=%f\n",iter,step_size);
/***********************************************************************************/


		int total_wirelength_verify = 0;
		int max_util = 0;
		int total_util = 0;
		int total_capacity = 0;
/*		printf("step size: %f\n", step_size);*/
		for (const auto &e : grid.edges()) {
			GridGraphEdgeObject &obj = get(edge_object, grid, e);
/*			assert(obj.max_util == 1);*/
			//obj.mult = std::max(0.0f, obj.mult+step_size*(obj.util-obj.max_util));
			//obj.slp=std::max(0.0f, (float) (obj.util-obj.max_util));
			//printf("slop or gradient : %f\n", obj.slp);
			obj.mult = obj.mult+step_size*(std::max(0.0f, (float) (obj.util-obj.max_util)));;
			//printf("util: %d mult: %f\n", obj.util, obj.mult);
			total_wirelength_verify += obj.util;
			max_util = std::max(max_util, obj.util);
			total_util += obj.util;
		}
/*************** nets size ********************************/
		//printf("Nets size is : %d\n", nets.size());
/**********************************************************/

		assert(total_wirelength_verify == total_wirelength);
		total_capacity = max_util * num_edges(grid);

		const float wire_delay = 0.5*101*2.874e-14 + 551*2.874e-14 + 58e-12;
//dhaarna-start
// putting values in the unordered maps enets and netse 
		/*timing analysis*/
		for (int i = 0; i < nets.size(); i++) {
/*			printf("Net: %s\n", nets[i].name.c_str());*/
/*			for (const auto &p : nets[i].points) {*/
/*				printf("%d,%d ", p.x, p.y);*/
/*			}*/
/*			printf("\n");*/

			Graph<adjacency_list<vecS, vecS, undirectedS, property<vertex_t_arr_t, float>, property<edge_weight_t, float>>, Point> route_tree;
			for (const auto &e : route_edges[i]) {Edge key;
Edges keys;
auto new_e=route_tree.add_edge(grid.object(source(e, grid)), grid.object(target(e, grid)));
if((grid.object(source(e, grid)).x<grid.object(target(e, grid)).x)||((grid.object(source(e, grid)).y<=grid.object(target(e, grid)).y) && (grid.object(source(e, grid)).x==grid.object(target(e, grid)).x)))
{Edge p(grid.object(source(e, grid)).x,grid.object(source(e, grid)).y, grid.object(target(e, grid)).x,grid.object(target(e, grid)).y);
Edges t(i,grid.object(source(e, grid)).x,grid.object(source(e, grid)).y, grid.object(target(e, grid)).x,grid.object(target(e, grid)).y); 
key=p;
keys=t;
				new_e = route_tree.add_edge(grid.object(source(e, grid)), grid.object(target(e, grid)));  }//dhaarna
else
{Edge p(grid.object(target(e, grid)).x,grid.object(target(e, grid)).y,grid.object(source(e, grid)).x,grid.object(source(e, grid)).y );  
Edges t(i, grid.object(target(e, grid)).x,grid.object(target(e, grid)).y,grid.object(source(e, grid)).x,grid.object(source(e, grid)).y); 
key=p;
keys=t;
				new_e = route_tree.add_edge(grid.object(target(e, grid)),grid.object(source(e, grid))); }//dhaarna
//std::cout<<grid.object(source(e, grid)).x<<grid.object(source(e, grid)).y<<grid.object(target(e, grid)).x<<grid.object(target(e, grid)).y<<" * "<<i<<":\n";
   enets[key].insert(i);  //enets stores nets corrs. to an edge
  netse[keys]=1;   //netse has value 1 if an edge is present in a net 


			put(edge_weight, route_tree, new_e, wire_delay);
			}
//dhaarna-end
			assert(is_tree(route_tree));
/*			printf("net: %d route_tree: %d\n", nets[i].points.size(), num_vertices(route_tree));*/
/*			assert(nets[i].points.size() == num_vertices(route_tree));*/
			vector<bool> visited(num_vertices(route_tree), false);
			get_delay(route_tree, route_tree.vertex(nets[i].points[0]), 0, visited);
/*		*/
			auto detailed_net = netlist.nets.find(nets[i].name)->second.get();
			assert(detailed_net->name == nets[i].name);
			
			for (const auto sink : detailed_net->sinks) {
				auto e = edge(detailed_net->source->tnode, sink->tnode, netlist.tgraph);
				assert(e.second == true);
/*				printf("%d,%d -> %d,%d: %f\n", detailed_net->source->port->block->position.x, detailed_net->source->port->block->position.y, */
/*						sink->port->block->position.x, sink->port->block->position.y,*/
/*						get(vertex_t_arr, route_tree, route_tree.vertex(sink->port->block->position)));*/
				put(edge_weight, netlist.tgraph, e.first, get(vertex_t_arr, route_tree, route_tree.vertex(sink->port->block->position)));
			}
		}


//dhaarna-start
//initial counts of edges violating 
int count1=0;
for (const auto &e : grid.edges()) {  

Edge key (grid.object(source(e, grid)).x,grid.object(source(e, grid)).y, grid.object(target(e, grid)).x,grid.object(target(e, grid)).y);                //dhaarna
if(enets[key].size()>MAX_UTIL){count1++;}
}
//std::cout<<"number of edges violating the constraint : "<<count1<<endl;

for (const auto &e : grid.edges()) {   
int np=0;        
Edge key;      
if((grid.object(source(e, grid)).x<grid.object(target(e, grid)).x)||((grid.object(source(e, grid)).y<=grid.object(target(e, grid)).y) && (grid.object(source(e, grid)).x==grid.object(target(e, grid)).x)))
{Edge p(grid.object(source(e, grid)).x,grid.object(source(e, grid)).y, grid.object(target(e, grid)).x,grid.object(target(e, grid)).y);
key=p;}
else
{Edge p(grid.object(target(e, grid)).x,grid.object(target(e, grid)).y,grid.object(source(e, grid)).x,grid.object(source(e, grid)).y );
key=p; }
//checking constraints
if(enets[key].size()>MAX_UTIL)
{
//std::cout<<"!!!!!!!!!!!!!!/n"<<enets[key].size()<<endl;
int d=enets[key].size() - MAX_UTIL;
int ax=grid.object(source(e, grid)).x;
int ay=grid.object(source(e, grid)).y;
int bx=grid.object(target(e, grid)).x;
int by=grid.object(target(e, grid)).y;
int tx,ty;
if(ax<bx || ((ax==bx) && ay<=by)){}
else {tx=ax;ty=ay;ax=bx;ay=by;bx=tx;by=ty;}

while(enets[key].size()>MAX_UTIL)
{
//bfs
//cout<<"#"<<enets[key].size()<<endl;
queue < pair<int,int> > q;
VisMap vis;
ParMap par;
int f=0;
pair<int,int> temp;
q.push(make_pair(ax,ay));
int iter=0;
while(!q.empty() && iter<=nx*ny)
{iter++;
temp=q.top();
if(temp==make_pair(bx,by)){f=1;break;}
q.pop();
for(int i=0;i<vic[temp].size();i++)
{if(vis[vic[temp][i]]!=1 )
{
Edge keyt;
if(temp.first<vic[temp][i].first|| ((temp.first==vic[temp][i].first) && temp.second<=vic[temp][i].second))
{Edge p(temp.first,temp.second,vic[temp][i].first,vic[temp][i].second); keyt=p;}
else
{Edge p(vic[temp][i].first,vic[temp][i].second,temp.first,temp.second); 
keyt=p;}
if(enets[keyt].size()<MAX_UTIL) 
{vis[vic[temp][i]]=1;
par[vic[temp][i]]=temp;
q.push(vic[temp][i]);

}
}}
}

Edge key(ax,ay,bx,by);

if(f==0)break;

else //path found
{
int thresh=INT_MAX;
pair<int,int> t=make_pair(bx,by);
while(t!=make_pair(ax,ay))
{
Edge w;
if(par[t].first<t.first|| ((par[t].first==t.first) && par[t].second<=t.second))
{Edge a(par[t].first,par[t].second,t.first,t.second); w=a;}
else
{
Edge a(t.first,t.second,par[t].first,par[t].second);w=a;
}
if(thresh> MAX_UTIL-enets[w].size())thresh=MAX_UTIL-enets[w].size();  // calculating thresh
t=par[t];
if(get<0>(key)== get<0>(w) && get<1>(key)==get<1>(w) && get<2>(key)==get<2>(w) && get<3>(key)==get<3>(w) ){break;}
}


vector <pair<int,int> > path;
 t=make_pair(bx,by);
path.push_back(t);
// forming path
while(t!=make_pair(ax,ay))
{t=par[t];
path.push_back(t);
}

//for(int r=0;r<path.size();r++)
//{
//std::cout<<ax <<ay<<bx<<by<<"pp : "<<path[r].first<<" "<<path[r].second<<endl;
//}


std::vector < pair<int,int> > v; // stores nets and the corrs count of edges that needs to be added to replace edge by path
for(auto e : enets[key])
{
int nn=e;
//replacing and updating data structures accordingly
int count=0;
for(int k=0;k<path.size()-1;k++)
{
Edge w;
Edges q;
if(path[k+1].first<path[k].first|| ((path[k+1].first==path[k].first) && path[k+1].second<=path[k].second))
{Edge a(path[k+1].first,path[k+1].second,path[k].first,path[k].second); w=a;
Edges b(nn,path[k+1].first,path[k+1].second,path[k].first,path[k].second); q=b;}
else
{
Edge a(path[k].first,path[k].second,path[k+1].first,path[k+1].second);w=a;
Edges b(nn,path[k].first,path[k].second,path[k+1].first,path[k+1].second);q=b;
}
if(get<0>(key)== get<0>(w) && get<1>(key)==get<1>(w) && get<2>(key)==get<2>(w) && get<3>(key)==get<3>(w) ){np=1;break;}
if(netse[q]!=1)
{                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
count++;
}
}
v.push_back(make_pair(count,nn));
}

std::sort(v.begin(),v.end());

int thresit;
for(thresit=0;thresit<thresh;thresit++ )
{
int nn=v[thresit].second;
//std::cout<<"nn : "<<nn<<endl;
list<GridGraph::edge_descriptor> :: iterator itrs;
itrs=route_edges[nn].begin();

while (itrs!= route_edges[nn].end()) {
int a1=grid.object(source(*itrs, grid)).x;
int b1=grid.object(source(*itrs, grid)).y;
int c1=grid.object(target(*itrs, grid)).x;
int d1=grid.object(target(*itrs, grid)).y;

if((ax==a1 && ay==b1 && bx==c1 && d1==by)||(ax==c1 && ay==d1 && bx==a1 && b1==by))
 {
route_edges[nn].erase(itrs);
break;
}

itrs++;}
//std::cout<<route_edges[nn].size()<<endl;
//std::cout<<enets[key].size()<<" "<<endl;
enets[key].erase(nn);
Edges keyst (nn,ax,ay,bx,by); //dhaarna
netse[keyst]=0;
//std::cout<<enets[key].size()<<" "<<endl;
//std::cout<<path.size()<<" "<<endl;
for(int k=0;k<path.size()-1;k++)
{
Edge w;
Edges q;
if(path[k+1].first<path[k].first|| ((path[k+1].first==path[k].first) && path[k+1].second<=path[k].second))
{Edge a(path[k+1].first,path[k+1].second,path[k].first,path[k].second); w=a;
Edges b(nn,path[k+1].first,path[k+1].second,path[k].first,path[k].second); q=b;}
else
{
Edge a(path[k].first,path[k].second,path[k+1].first,path[k+1].second);w=a;
Edges b(nn,path[k].first,path[k].second,path[k+1].first,path[k+1].second);q=b;
}
if(get<0>(key)== get<0>(w) && get<1>(key)==get<1>(w) && get<2>(key)==get<2>(w) && get<3>(key)==get<3>(w) ){np=1;break;}
if(netse[q]!=1)
{                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
//GridGraph g;

enets[w].insert(nn);

netse[q]=1;
GridGraph g;
GridGraph::edge_descriptor ep;
GridGraph::vertex_descriptor v1,v2;
v1=g.add_vertex(Point(path[k+1].first,path[k+1].second));
v2=g.add_vertex(Point(path[k].first,path[k].second));
if(path[k+1].first<path[k].first|| ((path[k+1].first==path[k].first) && path[k+1].second<=path[k].second))
{ep = g.add_edge(v1,v2);}
else
{ep = g.add_edge(v2,v1);}

route_edges[nn].push_back(ep);

}
}
d--;if(d==0 || np==1 ) break;}}if(np==1)break;}}

}
//wirelength verification
int wirel=0,wirel_verify=0;
for (int i = 0; i < nets.size(); i++) 
{wirel+=route_edges[i].size();}
int count=0;
for (const auto &e : grid.edges()) {  
//final count of edges violating 
Edge key (grid.object(source(e, grid)).x,grid.object(source(e, grid)).y, grid.object(target(e, grid)).x,grid.object(target(e, grid)).y);                //dhaarna
if(enets[key].size()>MAX_UTIL){count++;} 
if(enets[key].size() >= maxew)maxew=enets[key].size();
wirel_verify+=enets[key].size();
}



//if(wirel!=wirel_verify)std::cout<<"error in wirelength calc ! "<<endl;

//std::cout<<"number of edges violating the constraint : "<<count<<endl;
		//printf("---------\n");

		//while (current) {
			//printf("%s %s %d,%d %s[%d] arrival: %g\n", current->pin->port->block->name.c_str(), current->pin->port->block->instance.c_str(), current->pin->port->block->position.x, current->pin->port->block->position.y, current->pin->port->name.c_str(), current->pin->index, current->arrival_time);
			//current = current->prev;
		//}
		//printf("---------\n");
//execution time calculation 
		timer.stop();
		clock_t iter_end = clock();
		boost::timer::cpu_times elapsed = timer.elapsed();
		boost::timer::nanosecond_type time = elapsed.user + elapsed.system;
		boost::timer::nanosecond_type real_time = elapsed.wall;
// critical path delay calculation
		float average;
		timing_node crit_tnode;
		float crit = get_critical_path_delay(netlist.tgraph, &average, &crit_tnode);
		timing_node *current = &crit_tnode;
assert(wirel==wirel_verify);
//		printf("Iteration: %d Max_util: %d Total_wirelength: %d Runtime: %g, Crit_clb: %g Crit_T: %g Average_clb: %g Average_T: %g\n", iter, max_util, wirel_verify, (double)real_time/1000000000, crit, crit, average, average);
//std::cout<<"max ew for this iter "<<maxew<<endl;
total_routing_time += real_time;
		if (maxew <= minw) {
			minw=maxew;
			best_iter = iter;
			best_wirelength = wirel_verify;
			best_crit = crit;
			best_average = average;
		}
		if (crit < ult_crit) {
			ult_crit = crit;
			ult_crit_iter = iter;}}
	//printf("Total_runtime: %g Raw_time: %g Best_iteration: %d Util: %d Wirelength: %d Crit_clb: %g Crit_T: %g Average_clb: %g Average_T: %g\n", (double)total_routing_time/1000000000, (double)total_routing_time, best_iter, best_util, best_wirelength, best_crit, best_crit, best_average, best_average);
	//printf("Best crit: %g Best crit iteration: %d\n", ult_crit, ult_crit_iter);
	printf("%g,%d,%d,%g",(double)total_routing_time/1000000000,best_wirelength,minw ,ult_crit);
//std::cout<<"min w required : "<<minw<<endl;
}
//dhaarna- end
