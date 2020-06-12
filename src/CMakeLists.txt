# code relating to tree structured stick breaking process
include_directories(/usr/local/include/ ${CMAKE_SOURCE_DIR}/lib/eigen /usr/local/include/libhmsbeagle-1)

file( GLOB SRCS *.cpp *.hpp *.h)
add_library(sbn STATIC ${SRCS}) # create static library

find_package (PkgConfig REQUIRED)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

target_link_libraries(sbn ${PYTHON_LIBRARIES} /usr/local/lib/libhmsbeagle.dylib)
install(TARGETS sbn DESTINATION bin)


