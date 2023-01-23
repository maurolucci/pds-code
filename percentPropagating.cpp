//
// Created by max on 23.01.23.
//
#include <boost/program_options.hpp>
#include <iostream>
#include <fmt/format.h>
#include <range/v3/all.hpp>

#include "graphio.hpp"

int main(int argc, char** argv) {
    namespace po = boost::program_options;
    using std::string, std::vector;
    po::options_description desc;
    desc.add_options()
            ("graph,f", po::value<string>()->required(), "graph input file")
            ("outgraph,o", po::value<string>()->required(), "output graph")
            ("percentage,p", po::value<int>()->default_value(50), "percentage of propagating vertices to keep")
            ("all-non-prop,a", "consider all input vertices non-propagating")
            ("help,h", "show this message")
    ;
    po::positional_options_description pos;
    pos.add("graph", 1).add("outgraph", 1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cerr << "selects a subset of non-propagating vertices based on the non-propagating vertices in the input." << desc << std::endl;
        return 1;
    }
    auto graph = pds::readAutoGraph(vm["graph"].as<string>());
    if (vm.count("all-non-prop")) {
        for (auto v: graph.vertices()) {
            graph.getVertex(v).zero_injection = false;
        }
    }
    vector<pds::PowerGrid::VertexDescriptor> vertices;
    vertices.reserve(graph.numVertices());
    for (auto v: graph.vertices()) { if (!graph.getVertex(v).zero_injection) vertices.push_back(v); }
    ranges::shuffle(vertices);
    size_t numKeep = vertices.size() * vm["percentage"].as<int>() / 100;
    fmt::print("{} → {}\n", vm["graph"].as<string>(), vm["outgraph"].as<string>());
    fmt::print("keeping {} of {} non-propagating vertices", numKeep, vertices.size());
    for (size_t i = numKeep; i < vertices.size(); ++i) {
        graph.getVertex(vertices[i]).zero_injection = true;
    }
    size_t remaining = 0;
    for (auto v: graph.vertices()) { remaining += !graph.getVertex(v).zero_injection; }
    fmt::print(" ({})\n", remaining);
    pds::writePds(graph, vm["outgraph"].as<string>());
}