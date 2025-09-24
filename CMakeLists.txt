cmake_minimum_required(VERSION 2.8)

project(Router)

link_directories(
	/usr/boost_1_61_0/stage/lib	
	/usr/local/lib
	$ENV{TBB_LIB_DIR}
	)

set(srcs
	arch.cpp
	main.cpp
	netlist.cpp
	placement.cpp
	pugixml.cpp
	router.cpp
	)

set_source_files_properties(${srcs} PROPERTIES LANGUAGE CXX)

add_executable(Router ${srcs})

target_compile_options(Router PRIVATE
	-std=c++11
	-g -O3
	-Wall
	)

set(include_dirs
	/usr/boost_1_61_0
	/usr/local/include
	$ENV{TBB_INC_DIR}
	)

target_include_directories(Router PRIVATE
	${include_dirs}
	)

set(link_libs
	pthread
	dl 
	m
	boost_regex
	boost_timer
	boost_chrono
	boost_system
	tbb
	)

target_link_libraries(Router 
	${link_libs}
	)
