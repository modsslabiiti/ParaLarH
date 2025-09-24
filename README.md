Requirements:

C++ compiler support the C++11 standard
Boost C++ library -- install it using any package manager (eg. apt-get install libboost-all-dev) or build it from source (www.boost.org)
Intel TBB library -- download and extract the correct version for your platform from https://www.threadingbuildingblocks.org/download#stable-releases

Note*****Boost version 1.6***********
------------------------------------

Compilation steps:

1. cd build
2. export TBB_INC_DIR=<path to Intel TBB include directory>
   -- eg. export TBB_INC_DIR=/home/user/tbb/include (assuming the Intel TBB library has been extracted to /home/user/tbb)
3. export TBB_LIB_DIR=<path to Intel TBB library directory> 
   //-- eg. export TBB_INC_DIR=/home/user/tbb/lib (assuming the Intel TBB library has been extracted to /home/user/tbb)
	//-- eg. export TBB_LIB_DIR=/home/user/tbb/lib (assuming the Intel TBB library has been extracted to /home/user/tbb)
	this is working	   export TBB_LIB_DIR=/home/user/tbb/lib/intel64/gcc4.7
4. cmake ..			// generate make file folders
5. make				//Could you try running "make VERBOSE=1" and showing me the output?// gererate "Router executable file"
				// make clean command to clean make

After following the steps above, the final executable (Router) will be placed in the 'build' subdirectory.

------------------------------------

Before running the router:

The netlist required by our router has to be compiled using VTR (http://www.eecg.utoronto.ca/vtr/terms.html).

The command to do so is "./vpr <arch.xml> <benchmark_name> --pack --place".

//move to examples folder which have "vpr" executable file and right click on it tik mark it run as executable and from "examples" folder run:

eg. "./vpr fpga.xml ex5p --pack --place" would pack and place the benchmark called "ex5p" (in the file ex5p.blif) into the FPGA architecture described by the file "fpga.xml". The command would generate the files "ex5p.net" and "ex5p.place" which is required by our router. Note that "fpga.xml" and "ex5p.blif" has to be in the same directory as the "vpr" executable.
** VPR Default use Timing_Driven Router***
./vpr fpga.xml ex5p --pack --place --route
** to use Routing Driven Router **
./vpr fpga.xml ex5p --pack --place --route --router_algorithm breadth_first

In this package, an example FPGA architecture (fpga.xml) and a packed & placed benchmark (ex5p) have been included in the "examples" directory. Other benchmarks from the MCNC suite are also included but only in blif format which need to be compiled using VTR first.

------------------------------------
// Router in build directory....../home/rohita/paralar/build
//go in that location copy it and paste it to examples directory
Running the router:
./Router <arch.xml> <benchmark_name> <num_threads> <num_of_iterations> <max_channel_capacity>

eg. "./Router fpga.xml ex5p 4 50 100" would route the "ex5p" benchmark with a channel width of 100 using 4 threads for 50 iterations. Note that the files "fpga.xml", "ex5p.blif", "ex5p.net", and "ex5p.place" has to be in the same directory as the "Router" executable.




gdb --args ./Router <arch.xml> <benchmark_name> <num_threads> <num_of_iterations> <max_channel_capacity>


