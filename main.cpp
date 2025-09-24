#include "pch.hpp"
#include "mst.hpp"
#include "netlist.hpp"
#include "router.hpp"
#include "arch.hpp"
#include <boost/version.hpp>

int MAX_ITER;
int MAX_UTIL;

int main(int argc, char *argv[])
{
	
printf("boost ver %d\n",BOOST_VERSION);
assert(argc == 6);

	int num_threads = atoi(argv[3]);
	MAX_ITER = atoi(argv[4]);
	MAX_UTIL = atoi(argv[5]);

	tbb::task_scheduler_init init(num_threads);

	Arch arch;
	arch.load(argv[1]);

	Netlist n;
	n.load(argv[2], arch.pbs);

	Router router(n.grid_size.first, n.grid_size.second);

	router.parallel_route_tbb_new_netlist(n, num_threads);

	return 0;
}
