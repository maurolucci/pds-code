#include <iostream>

#include <boost/graph/graphml.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/property_map/compose_property_map.hpp>

#include <string>

#include "pds.hpp"
#include "gurobi_solve.hpp"

template<class T, class Graph>
auto make_vector_property_map(T&& vec, Graph graph) {
    return boost::make_iterator_property_map(vec.begin(), get(boost::vertex_index, graph));
}

template<class ActiveMap, class ObservedMap>
void printResult(const pds::PowerGrid& graph, const ActiveMap& active, const ObservedMap& observed) {
    bool sep = false;
    size_t active_count = 0, inactive_count = 0, zero_injection_count = 0;
    std::cout << "active nodes: [";
    for (const auto &v: graph.vertices()) {
        switch (active[v]) {
            case pds::PmuState::Active:
                if (sep) std::cout << ", ";
                sep = true;
                std::cout << graph[v].name;
                ++active_count;
                break;
            case pds::PmuState::Inactive:
                ++inactive_count;
                break;
            case pds::PmuState::Blank:
                break;
        }
        zero_injection_count += graph[v].zero_injection;
    }
    std::cout << "]\n";
    std::cout << "graph (m=" << graph.num_edges() << ", n=" << graph.num_vertices()
              << ", #active=" << active_count
              << ", #inactive=" << inactive_count
              << ", #zero_injection=" << zero_injection_count << ")\n";
    std::cout << "feasible: " << std::boolalpha << pds::observed(graph, observed) << std::endl;
}

void processBoost(const std::string& filename) {
    auto graph = pds::import_graphml(filename);
    for (auto it = vertices(graph).first; it != vertices(graph).second; ++it) {
        auto v = *it;
        graph[v].zero_injection = true;//zero_injection[v];
    }
    {
        pds::map<decltype(graph)::vertex_descriptor, pds::PmuState> active_data;
        auto active = boost::associative_property_map(active_data);
        pds::solve_pds(graph, active);
        pds::map<decltype(graph)::vertex_descriptor, bool> observed_data;
        pds::dominate(graph, active, boost::associative_property_map(observed_data));
        pds::propagate(graph, graph.getProperty(&pds::Bus::zero_injection), boost::associative_property_map(observed_data));
        printResult(graph, active, boost::associative_property_map(observed_data));
    }
    {
        pds::PowerGrid grid(graph);
        pds::map<std::string, pds::PmuState> active_data;
        {
            auto active = boost::compose_property_map(boost::associative_property_map(active_data), get(&pds::Bus::name, grid));
            pds::map<std::string, bool> observed_data;
            auto observed = boost::compose_property_map(boost::associative_property_map(observed_data), get(&pds::Bus::name, grid));
            auto zero_injection = get(&pds::Bus::zero_injection, graph);
            std::cout << "before presolving: " << std::endl;
            printResult(grid, active, observed);
            pds::presolve(grid, zero_injection, active, observed);
            std::cout << "after presolving: " << std::endl;
            printResult(grid, active, observed);
        }
        auto active = boost::compose_property_map(boost::associative_property_map(active_data), get(&pds::Bus::name, graph));
        pds::solve_pds(grid, active);
        pds::map<std::string, bool> observed_map;
        auto observed = boost::compose_property_map(boost::associative_property_map(observed_map), get(&pds::Bus::name, graph));
        pds::dominate(graph, active, observed);
        pds::propagate(graph, get(&pds::Bus::zero_injection, graph),observed);
        printResult(graph, active, observed);
    }
    if (false) {
        std::cout << "take 2:" << std::endl;
        pds::map<decltype(graph)::vertex_descriptor, pds::PmuState> active_data;
        auto active = boost::associative_property_map(active_data);
        pds::map<decltype(graph)::vertex_descriptor, bool> observed_data;
        auto observed = boost::associative_property_map(observed_data);
        for (const auto& v: pds::range_pair(vertices(graph))) {
            std::vector<int> pmus{5, 11, 55};
            if (ranges::contains(
                    pmus | ranges::views::transform([](int i) -> std::string { return std::to_string(i);}),
                    graph[v].name)
            ) {
                active[v] = pds::PmuState::Active;
            }
        }
        pds::dominate(graph, active, observed);
        pds::propagate(graph, get(&pds::Bus::zero_injection, graph), observed);
        printResult(graph, active, observed);
    }
}

int main(int argc, const char** argv) {
    std::string filename = argv[1];
    processBoost(filename);
    return 0;
}
