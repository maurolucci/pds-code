cmake_minimum_required(VERSION 3.22)
project(pds)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 20)

find_package(Boost 1.70 REQUIRED COMPONENTS program_options)

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
FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG c4ee726 # 9.0.0
)
FetchContent_Declare(
        htd
        GIT_REPOSITORY https://github.com/CinusMinus/htd.git
        GIT_TAG master
)
FetchContent_MakeAvailable(unordered_dense OGDF fmt htd)
# No warnings from OGDF headers
get_target_property(TMP_OGDF_INCLUDE_DIRS OGDF INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(OGDF SYSTEM INTERFACE ${TMP_OGDF_INCLUDE_DIRS})

add_subdirectory(contrib/mpgraphs)
find_package(Gurobi REQUIRED)
find_package(range-v3 REQUIRED)
find_package(tinyxml2 REQUIRED)

function(enable_warnings target)
    target_compile_options(${target} PRIVATE
            $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
            -Wall -Wextra -Werror>
            $<$<CXX_COMPILER_ID:MSVC>:
            /W4>)
endfunction()

add_library(pdslib STATIC src/pds/pds.cpp src/pds/graphml.cpp src/pds/draw_grid.cpp src/pds/graphio.cpp)
target_include_directories(pdslib PUBLIC include)
target_link_libraries(pdslib PRIVATE tinyxml2)
target_link_libraries(pdslib PUBLIC range-v3 fmt::fmt)
target_link_libraries(pdslib PUBLIC mpgraphs unordered_dense::unordered_dense PRIVATE OGDF)
enable_warnings(pdslib)

add_library(pdslib-gurobi src/pds/gurobi_solve.cpp)
target_link_libraries(pdslib-gurobi PUBLIC pdslib PRIVATE gurobi::gurobi_cxx)
enable_warnings(pdslib-gurobi)

add_executable(pds main.cpp)
target_link_libraries(pds PUBLIC pdslib pdslib-gurobi fmt::fmt Boost::program_options)
enable_warnings(pds)

add_executable(random_tree random_trees.cpp)
target_link_libraries(random_tree PUBLIC pdslib Boost::program_options)

add_executable(shave shave_graph.cpp)
target_link_libraries(shave PUBLIC pdslib fmt::fmt Boost::program_options)

add_executable(check_graph check_graph.cpp)
target_link_libraries(check_graph PUBLIC pdslib fmt::fmt Boost::program_options)

add_executable(experiment experiment.cpp)
target_link_libraries(experiment PUBLIC pdslib pdslib-gurobi fmt::fmt Boost::program_options htd)

add_executable(draw_grid drawGrid.cpp)
target_link_libraries(draw_grid PUBLIC pdslib fmt::fmt Boost::program_options OGDF)

add_subdirectory(testgraphs)