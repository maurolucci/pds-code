cmake_minimum_required(VERSION 3.22)
project(pds)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
FetchContent_Declare(
        ogdf
        URL https://ogdf.uos.de/wp-content/uploads/2022/02/ogdf.v2022.02.zip
)
FetchContent_MakeAvailable(ogdf)

find_package(Boost 1.70 REQUIRED COMPONENTS graph)

add_library(pdslib INTERFACE)
target_include_directories(pdslib INTERFACE include)

find_package(Gurobi REQUIRED)
find_package(range-v3 REQUIRED)

add_executable(pds main.cpp include/pds.hpp include/gurobi_solve.cpp include/gurobi_solve.hpp)

add_library(graphml src/graphml.cpp)

target_link_libraries(pds PUBLIC OGDF)
target_link_libraries(pds PUBLIC Boost::graph range-v3 gurobi::gurobi_cxx)
target_link_libraries(pds PUBLIC pdslib)
