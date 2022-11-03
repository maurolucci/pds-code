//
// Created by max on 15.08.22.
//

#include <mpgraphs/graph.hpp>
#include <mpgraphs/boost_adapter.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <range/v3/all.hpp>
#include <boost/graph/graphml.hpp>

#include <boost/program_options.hpp>

#include <fmt/format.h>

#include "map.hpp"

using Graph = mpgraphs::MapGraph<mpgraphs::Empty, mpgraphs::Empty, mpgraphs::EdgeDirection::Undirected>;
using Vertex = Graph::VertexDescriptor;

template<std::uniform_random_bit_generator Rng = std::minstd_rand>
Graph randomTree(size_t n, Rng rng = Rng()) {
    Graph graph;
    if (n == 0) { return graph; }
    std::vector<Vertex> graphVertices;
    graphVertices.push_back(graph.addVertex());
    for (size_t x = 1; x < n; ++x) {
        auto u = graphVertices[std::uniform_int_distribution(size_t{0}, graphVertices.size() - 1)(rng)];
        auto v = graph.template addVertex();
        graphVertices.push_back(v);
        graph.addEdge(u, v);
    }
    return graph;
}

template<typename L, typename R>
bool intersects(const L& lhs, const R& rhs) {
    if (lhs.size() < rhs.size()) {
        intersects(rhs, lhs);
    } else {
        for (auto x: rhs) {
            if (lhs.contains(x)) return true;
        }
    }
    return false;
}

void dfs(const Graph& graph, Vertex start, pds::map<Vertex, Vertex>& parent, pds::map<Vertex, size_t>& distance) {
    if (!parent.contains(start)) {
        parent[start] = start;
        distance[start] = 0;
    }
    for (auto w: graph.neighbors(start)) {
        if (!parent.contains(w)) {
            parent[w] = start;
            distance[w] = distance[start] + 1;
            dfs(graph, w, parent, distance);
        }
    }
}

std::optional<Vertex> lowestCommonAncestor(Vertex first, Vertex second, const pds::map<Vertex, Vertex>& parent, const pds::map<Vertex, size_t>& distance) {
    while (first != second) {
        if (parent.at(first) == first && parent.at(second) == second) return {};
        if ((distance.at(first) < distance.at(second) && parent.at(second) != second) || parent.at(first) == first) {
            second = parent.at(second);
        } else {
            first = parent.at(first);
        }
    }
    return first;
}

template<std::uniform_random_bit_generator Rng = std::minstd_rand>
Graph randomCactus(size_t n, Rng rng = Rng()) {
    Graph graph;
    if (n == 0) { return graph; }
    std::vector<Vertex> graphVertices;
    std::set<Vertex> blocks;
    graphVertices.push_back(graph.addVertex());
    n -= 1;
    while (n > 0)  {
        auto u = graphVertices[std::uniform_int_distribution(size_t{0}, graphVertices.size() - 1)(rng)];
        auto v = graph.template addVertex();
        graphVertices.push_back(v);
        graph.addEdge(u, v);
        if (!blocks.contains(u)) {
            blocks.insert(v);
        } else {
            --n;
        }
    }
    for (auto v: blocks) {
        auto neighbors = graph.neighbors(v) | ranges::to<std::vector>();
        graph.removeVertex(v);
        ranges::shuffle(neighbors, rng);
        auto d = neighbors.size();
        for (size_t i = 0; i < neighbors.size(); ++i) {
            graph.addEdge(neighbors[i], neighbors[(i+1) % d]);
        }
    }
    return graph;
}

//Graph randomCactus(size_t n, double p = 0.1) {
//    auto graph = randomTree(n);
//    if (n == 0) return graph;
//    pds::map<Graph::vertex_descriptor, pds::set<size_t>> onCycle;
//    std::minstd_rand rng;
//    const auto N = graph.numVertices();
//    std::uniform_int_distribution endpoint(size_t{0}, N - 1);
//    auto vertices = graph.vertices() | ranges::to<std::vector>();
//    pds::map<Vertex, Vertex> parent;
//    pds::map<Vertex, size_t> distance;
//    dfs(graph, vertices.front(), parent, distance);
//    size_t cycle = 1;
//    for (size_t i = 0; i < static_cast<size_t>(N * p); ++i) {
//        auto u = endpoint(rng);
//        auto v = endpoint(rng);
//        auto x = vertices[u];
//        auto y = vertices[v];
//        if (!intersects(onCycle[x], onCycle[y])) {
//            auto lca = lowestCommonAncestor(x, y, parent, distance).value();
//            graph.addEdge(x, y);
//            while (x != lca) {
//                onCycle[x].insert(cycle);
//                x = parent[x];
//            }
//            while (y != lca) {
//                onCycle[y].insert(cycle);
//                y = parent[y];
//            }
//            onCycle[lca].insert(cycle);
//            ++cycle;
//        }
//    }
//    return graph;
//}

namespace po = boost::program_options;

void writeGraph(const Graph& graph, std::ostream& out) {
    pds::map<Vertex, long> zeroInjection;
    for (auto v: graph.vertices()) {
        zeroInjection[v] = 1;
    }
    boost::associative_property_map ziMap(zeroInjection);
    boost::dynamic_properties properties(boost::ignore_other_properties);
    properties.property("zero_injection", ziMap);
    pds::map<Vertex, size_t> vertexIndex;
    for (size_t i = 0; auto v: graph.vertices()) {
        vertexIndex.insert({v, i});
        ++i;
    }
    boost::write_graphml(out, graph, boost::associative_property_map(vertexIndex), properties);
}

int main(int argc, char** argv) {
    po::options_description desc("options");
    std::string graphType;
    std::string outFile;
    size_t n = 0;
    desc.add_options()
            ("help,h", "show this help")
            ("type,t", po::value(&graphType)->default_value("tree"), "graph type")
            ("size,n", po::value(&n)->default_value(100), "number of vertices")
            ("seed,s", po::value<size_t>(), "random seed")
            ("out,o", po::value(&outFile)->default_value("-"), "out file name")
            ;
    po::positional_options_description positional;
    positional.add("size", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(positional).run(), vm);
    po::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    Graph graph;
    std::minstd_rand rng(vm.count("seed") ? vm["seed"].as<size_t>() : size_t{std::random_device()()});
    if (graphType == "tree") { graph = randomTree(n, rng);}
    else if (graphType == "cactus") { graph = randomCactus(n, rng); }

    std::optional<std::ofstream> outStream;
    if (outFile != "-") {
        outStream.emplace(std::ofstream(outFile));
    }
    writeGraph(graph, outStream.has_value() ? outStream.value() : std::cout);
}