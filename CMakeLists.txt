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

function(bito_compile_opts PRODUCT)
	target_compile_features(${PRODUCT} PUBLIC cxx_std_17)
	target_compile_options(${PRODUCT} PUBLIC -Werror)
	target_include_directories(${PRODUCT} PUBLIC
		${PROJECT_BINARY_DIR}/beagle-lib/install/include/libhmsbeagle-1
		lib/eigen
	)
endfunction()

function(bito_link_opts PRODUCT)
	target_link_options(${PRODUCT} PUBLIC -Wl,-rpath,.)
	target_link_libraries(${PRODUCT} PUBLIC bito-core)
	add_dependencies(${PRODUCT} bito-core)
endfunction()

function(bito_executable PRODUCT)
    add_executable(${PRODUCT} ${ARGN})
	bito_compile_opts(${PRODUCT})
	bito_link_opts(${PRODUCT})
endfunction()

###################
# libbito-core.so #
###################

add_library(bito-core SHARED
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
	src/subsplit_dag_nni.cpp
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
bito_compile_opts(bito-core)
target_link_options(bito-core PUBLIC -pthread -Wl,-rpath,beagle-lib/install/lib)
target_link_libraries(bito-core PUBLIC hmsbeagle hmsbeagle-cpu z)
target_link_directories(bito-core PUBLIC
	${PROJECT_BINARY_DIR}/beagle-lib/install/lib
)
add_dependencies(bito-core beagle-lib)

##############
# libbito.so #
##############

add_library(bito SHARED
	src/pybito.cpp
)
bito_compile_opts(bito)
target_include_directories(bito PUBLIC
	${Python_INCLUDE_DIRS}
	${PROJECT_BINARY_DIR}/pybind11/install/include
)
bito_link_opts(bito)
add_dependencies(bito pybind11)

###############
# executables #
###############

bito_executable(doctest
	src/doctest.cpp)

bito_executable(gp_doctest
	src/gp_doctest.cpp)

bito_executable(noodle
	src/noodle.cpp)

set(CMAKE_SKIP_RPATH ON)
