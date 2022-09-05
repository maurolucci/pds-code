//
// Created by max on 18.08.22.
//

#include <setgraph/graph.hpp>
#include <range/v3/all.hpp>
#include <boost/program_options.hpp>

#include "pds.hpp"
#include "graphio.hpp"

using Graph = pds::PowerGrid;//setgraph::SetGraph<setgraph::Empty, setgraph::Empty, setgraph::EdgeDirection::Undirected>;
using Vertex = Graph::vertex_descriptor;
using Edge = Graph::edge_descriptor;

void shave(Graph& graph, bool leaves = true, bool paths = true) {
    std::vector<Vertex> deg1;
    std::vector<Vertex> deg2;
    for (Vertex v: graph.vertices()) {
        graph.removeEdge(v, v);
        if (graph.degree(v) <= 1) {
            deg1.push_back(v);
        } else if (graph.degree(v) == 2) {
            deg2.push_back(v);
        }
    }
    while ((leaves && !deg1.empty()) || (paths && !deg2.empty())) {
        while (leaves && !deg1.empty()) {
            auto v = deg1.back(); deg1.pop_back();
            if (graph.hasVertex(v)) {
                for (auto w: graph.neighbors(v)) {
                    if (graph.degree(w) <= 2) {
                        deg1.push_back(w);
                    } else if (graph.degree(w) == 3) {
                        deg2.push_back(w);
                    }
                }
                graph.removeVertex(v);
            }
        }
        while (paths && !deg2.empty()) {
            auto v = deg2.back(); deg2.pop_back();
            if (graph.hasVertex(v)) {
                if (graph.degree(v) == 2) {
                    auto neighbors = graph.neighbors(v) | ranges::to<std::vector>();
                    if (!graph.edge(neighbors[0], neighbors[1])) {
                        graph.addEdge(neighbors[0], neighbors[1]);
                        graph.removeVertex(v);
                    }
                } else if (graph.degree(v) < 2) {
                    deg1.push_back(v);
                }
            }
        }
    }
}

int main(int argc, const char** argv) {
    namespace po = boost::program_options;
    po::options_description desc;
    using std::string, std::vector;
    desc.add_options()
            ("graph,f", po::value<string>(), "graph input file")
            ("outfile,o", po::value<string>()->default_value("-"), "graph output file")
            ("leaves,l", "remove leaves")
            ("paths,p","contract paths")
            ;
    po::positional_options_description pos;
    pos.add("graph", 1).add("outfile", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
    po::notify(vm);

    auto graph = pds::readGraphML(vm["graph"].as<string>());
    shave(graph, vm.count("leaves"), vm.count("paths"));

    std::optional<std::ofstream> outStream;
    auto outfile = vm["outfile"].as<string>();
    if (outfile != "-") {
        outStream = {std::ofstream(outfile)};
    }
    pds::exportGraphml(graph, outStream.has_value() ? *outStream : std::cout);
}
