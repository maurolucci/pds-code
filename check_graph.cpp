//
// Created by max on 19.08.22.
//

#include <setgraph/graph.hpp>
#include <range/v3/all.hpp>
#include <boost/program_options.hpp>

#include "pds.hpp"
#include "graphio.hpp"

namespace pds {

using Vertex = PowerGrid::vertex_descriptor;

bool visitComponent(const PowerGrid& graph, Vertex start, set<Vertex>& seen) {
    if (seen.contains(start)) return false;
    seen.insert(start);
    for (auto w: graph.neighbors(start)) {
        visitComponent(graph, w, seen);
    }
    return true;
}

size_t numComponents(const PowerGrid &graph) {
    set<Vertex> seen;
    size_t nComponents = 0;
    for (auto v: graph.vertices()) {
        if (visitComponent(graph, v, seen)) ++nComponents;
    }
    return nComponents;
}

bool isConnected(const PowerGrid& graph) {
    set<Vertex> seen;
    size_t nComponents = 0;
    for (auto v: graph.vertices()) {
        if (visitComponent(graph, v, seen)) ++nComponents;
        if (nComponents > 1) return false;
    }
    return true;
}

bool isTree(const PowerGrid& graph) {
    return graph.numEdges() + 1 == graph.numVertices() && isConnected(graph);
}

bool isForest(const PowerGrid& graph) {
    return graph.numEdges() + numComponents(graph) == graph.numVertices();
}

}

int main(int argc, const char** argv) {
    using namespace pds;
    auto graph = pds::readGraphML(argv[1]);
    fmt::print("n={}, m={}\n", graph.numVertices(), graph.numEdges());
    fmt::print("is connected: {}\n", isConnected(graph));
    fmt::print("is forest: {}\n", isForest(graph));
    fmt::print("is tree: {}\n", isTree(graph));
}
