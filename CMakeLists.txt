cmake_minimum_required(VERSION 3.22)
project(pds)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 20)

find_package(Boost 1.70 REQUIRED COMPONENTS graph)

include(FetchContent)
set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE OPT_IN)
FetchContent_Declare(
        unordered_dense
        GIT_REPOSITORY https://github.com/martinus/unordered_dense.git
        GIT_TAG 94156450d991515faa07ccd3823a2592c000e49f
)
FetchContent_Declare(
        OGDF
        URL https://ogdf.uos.de/wp-content/uploads/2022/02/ogdf.v2022.02.zip
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(unordered_dense OGDF)

add_library(pdslib STATIC src/pds.cpp src/graphml.cpp src/draw_grid.cpp)
target_include_directories(pdslib PUBLIC include)
target_compile_options(pdslib PRIVATE
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
        -Wall -Wextra>
        $<$<CXX_COMPILER_ID:MSVC>:
        /W4>)

add_subdirectory(contrib/setgraph)
target_link_libraries(pdslib PUBLIC setgraph unordered_dense::unordered_dense PRIVATE OGDF)

find_package(Gurobi REQUIRED)
find_package(range-v3 REQUIRED)

add_executable(pds main.cpp include/pds.hpp include/gurobi_solve.cpp include/gurobi_solve.hpp include/map.hpp include/draw_grid.hpp)

target_link_libraries(pds PUBLIC Boost::graph range-v3 gurobi::gurobi_cxx)
target_link_libraries(pds PUBLIC pdslib)