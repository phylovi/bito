cmake_minimum_required(VERSION 3.19)
project(bito)
find_package(Python COMPONENTS Development)

include(ExternalProject)
ExternalProject_Add(beagle-lib
	GIT_REPOSITORY	https://github.com/beagle-dev/beagle-lib.git
	GIT_TAG			origin/hmc-clock
	GIT_SHALLOW		true
	GIT_PROGRESS	true
	CMAKE_ARGS		-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/beagle-lib/install -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF -DBUILD_JNI=OFF
	PREFIX			beagle-lib
	INSTALL_DIR		beagle-lib/install
)

include(ExternalProject)
ExternalProject_Add(pybind11
	GIT_REPOSITORY	https://github.com/pybind/pybind11.git
	GIT_TAG			origin/v2.7
	GIT_SHALLOW		true
	GIT_PROGRESS	true
	CMAKE_ARGS		-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/pybind11/install -DPYBIND11_TEST=OFF
	PREFIX			pybind11
	INSTALL_DIR		pybind11/install
)

add_library(bito SHARED
	src/alignment.cpp
	src/bitset.cpp
	src/block_model.cpp
	src/block_specification.cpp
	src/clock_model.cpp
	src/combinatorics.cpp
	src/csv.cpp
	src/driver.cpp
	src/engine.cpp
	src/fat_beagle.cpp
	src/gp_dag.cpp
	src/gp_engine.cpp
	src/gp_instance.cpp
	src/gp_operation.cpp
	src/mersenne_twister.cpp
	src/node.cpp
	src/numerical_utils.cpp
	src/parser.cpp
	src/phylo_model.cpp
	src/psp_indexer.cpp
	src/quartet_hybrid_request.cpp
	src/rooted_gradient_transforms.cpp
	src/rooted_sbn_instance.cpp
	src/rooted_tree.cpp
	src/rooted_tree_collection.cpp
	src/sbn_maps.cpp
	src/sbn_probability.cpp
	src/sbn_support.cpp
	src/scanner.cpp
	src/site_model.cpp
	src/site_pattern.cpp
	src/stick_breaking_transform.cpp
	src/subsplit_dag.cpp
	src/subsplit_dag_node.cpp
	src/substitution_model.cpp
	src/taxon_name_munging.cpp
	src/tidy_subsplit_dag.cpp
	src/tree.cpp
	src/tree_collection.cpp
	src/unrooted_sbn_instance.cpp
	src/unrooted_tree.cpp
	src/unrooted_tree_collection.cpp
	src/zlib_stream.cpp
)
target_compile_options(bito PUBLIC -std=c++17)
target_link_options(bito PUBLIC -pthread -Wl,-rpath,beagle-lib/install/lib)
target_link_libraries(bito PUBLIC hmsbeagle hmsbeagle-cpu z)
target_include_directories(bito PUBLIC
	lib/eigen
	${PROJECT_BINARY_DIR}/beagle-lib/install/include/libhmsbeagle-1
)
target_link_directories(bito PUBLIC
	${PROJECT_BINARY_DIR}/beagle-lib/install/lib
)
add_dependencies(bito beagle-lib)



add_library(bito-python SHARED
	src/pylibsbn.cpp
)
target_compile_options(bito-python PUBLIC -std=c++17)
target_include_directories(bito-python PUBLIC
	${Python_INCLUDE_DIRS}
	${PROJECT_BINARY_DIR}/pybind11/install/include
	${PROJECT_BINARY_DIR}/beagle-lib/install/include/libhmsbeagle-1
	lib/eigen
)
add_dependencies(bito-python bito pybind11)

set(CMAKE_SKIP_RPATH ON)
