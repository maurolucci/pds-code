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

struct Bus {
    bool zero_injection;
};

template<class T>
class WrapBool {
    bool& source;
public:
    explicit WrapBool(bool& source) : source(source) { }
    WrapBool& operator=(const T& other) {
        source = bool(other);
        return *this;
    }
    operator T() {
        return T(source);
    }
};

template<class T, class Graph>
auto make_vector_property_map(T&& vec, Graph graph) {
    return boost::make_iterator_property_map(vec.begin(), get(boost::vertex_index, graph));
}

void processBoost(const std::string& filename)
{
    std::ifstream graph_in(filename);
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Bus> graph;
    using vertex_id_t = decltype(graph)::vertex_descriptor;
    using edge_id_t = decltype(graph)::edge_descriptor;
    //std::vector<long> long_zi(boost::num_vertices(graph));
    boost::dynamic_properties attr(boost::ignore_other_properties);
    std::vector<long> long_zi(num_vertices(graph));
    typename boost::property_map<decltype(graph), boost::vertex_index_t>::type index = get(boost::vertex_index, graph);
    boost::vector_property_map<long, decltype(index)> zero_injection(num_vertices(graph), index);
    attr.property("zero_injection", zero_injection);//make_vector_property_map(long_zi, graph));
    boost::read_graphml(graph_in, graph, attr);
    graph_in.close();
    for (auto it = vertices(graph).first; it != vertices(graph).second; ++it) {
        auto v = *it;
        graph[v].zero_injection = true;//zero_injection[v];
    }
    std::cout << "n=" << boost::num_vertices(graph) << ", m=" << boost::num_edges(graph) << std::endl;
    std::cout << "attributes: " << std::distance(attr.begin(), attr.end()) << std::endl;
    for (auto a: attr) {
        std::cout << "found attribute " << a.first << std::endl;
    }
    //std::vector<pds::PmuState> active(boost::num_vertices(graph));
    //std::vector<uint8_t> observed(boost::num_vertices(graph));
    //auto active_map = make_vector_property_map(active, graph);
    //auto observed_map = make_vector_property_map(observed, graph);
    {
        boost::vector_property_map<pds::PmuState, decltype(index)> active_map(num_vertices(graph), index);
        boost::vector_property_map<uint8_t, decltype(index)> observed_map(num_vertices(graph), index);
        pds::solve_pds(graph, active_map, observed_map);
        pds::dominate(graph, active_map, observed_map);
        {
            bool sep = false;
            std::cout << "active vertices: ";
            for (const auto &v: pds::range_pair(vertices(graph))) {
                if (active_map[v] == pds::PmuState::Active) {
                    if (sep) std::cout << ", ";
                    sep = true;
                    std::cout << v;
                }
            }
        }
        std::cout << std::endl;
        pds::propagate(graph, observed_map);
        std::cout << "fully observed: " << pds::observed(graph, observed_map) << std::endl;
    }

    {
        std::cout << "take 2:" << std::endl;
        boost::vector_property_map<pds::PmuState, decltype(index)> active_map(num_vertices(graph), index);
        boost::vector_property_map<uint8_t, decltype(index)> observed_map(num_vertices(graph), index);
        for (auto i: std::vector{5, 11, 55}) {
            active_map[i] = pds::PmuState::Active;
        }
        pds::dominate(graph, active_map, observed_map);
        pds::propagate(graph, observed_map);
        std::cout << "fully observed: " << pds::observed(graph, observed_map) << std::endl;
    }
}

int main(int argc, const char** argv) {
    std::cout << "Hello, World!" << std::endl;
    std::string filename = argv[1];
    processOgdf(filename);
    processBoost(filename);
    return 0;
}
