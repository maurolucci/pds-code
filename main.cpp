#include <iostream>
#include "ogdf/basic/Graph.h"
#include "ogdf/fileformats/GraphIO.h"

#include <boost/graph/graphml.hpp>
#include <boost/property_map/function_property_map.hpp>

#include "pds.hpp"
#include "gurobi_solve.hpp"

void processOgdf(const std::string& filename) {
    ogdf::Graph graph;
    ogdf::GraphAttributes attr;
    std::cout << filename << std::endl;
    std::cout << "success: " << ogdf::GraphIO::read(attr, graph, filename) << std::endl;
    std::cout << graph.numberOfNodes() << " " << graph.numberOfEdges() << std::endl;
    auto first_node = graph.firstNode();
    std::cout << "attributes: " << attr.attributes() << std::endl;
}

template<class T, class Graph>
auto make_vector_property_map(T&& vec, Graph graph) {
    return boost::make_iterator_property_map(vec.begin(), get(boost::vertex_index, graph));
}

void processBoost(const std::string& filename)
{
    auto graph = pds::read_graphml(filename);
    for (auto it = vertices(graph).first; it != vertices(graph).second; ++it) {
        auto v = *it;
        graph[v].zero_injection = true;//zero_injection[v];
    }
    std::cout << "n=" << boost::num_vertices(graph) << ", m=" << boost::num_edges(graph) << std::endl;
    //std::vector<pds::PmuState> active(boost::num_vertices(graph));
    //std::vector<uint8_t> observed(boost::num_vertices(graph));
    //auto active_map = make_vector_property_map(active, graph);
    //auto observed_map = make_vector_property_map(observed, graph);
    {
        boost::associative_property_map<std::map<decltype(graph)::vertex_descriptor, pds::PmuState>> active;
        boost::associative_property_map<std::map<decltype(graph)::vertex_descriptor, bool>> observed;
        pds::solve_pds(graph, active, observed);
        pds::dominate(graph, active, observed);
        {
            bool sep = false;
            std::cout << "active vertices: ";
            for (const auto &v: pds::range_pair(vertices(graph))) {
                if (get(active, v) == pds::PmuState::Active) {
                    if (sep) std::cout << ", ";
                    sep = true;
                    std::cout << v;
                }
            }
        }
        std::cout << std::endl;
        pds::propagate(graph, observed);
        std::cout << "fully observed: " << pds::observed(graph, observed) << std::endl;
    }

    {
        std::cout << "take 2:" << std::endl;
        boost::associative_property_map<std::map<decltype(graph)::vertex_descriptor, pds::PmuState>> active;
        boost::associative_property_map<std::map<decltype(graph)::vertex_descriptor, bool>> observed;
        for (auto i: std::vector{5, 11, 55}) {
            active[i] = pds::PmuState::Active;
        }
        pds::dominate(graph, active, observed);
        pds::propagate(graph, observed);
        std::cout << "fully observed: " << pds::observed(graph, observed) << std::endl;
    }
}

int main(int argc, const char** argv) {
    std::cout << "Hello, World!" << std::endl;
    std::string filename = argv[1];
    processOgdf(filename);
    processBoost(filename);
    return 0;
}
