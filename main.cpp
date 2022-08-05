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
        pds::map<long, pds::PmuState> active_data;
        {
            auto active = boost::compose_property_map(boost::associative_property_map(active_data), get(&pds::Bus::id, grid));
            pds::map<long, bool> observed_data;
            auto observed = boost::compose_property_map(boost::associative_property_map(observed_data), get(&pds::Bus::id, grid));
            auto zero_injection = get(&pds::Bus::zero_injection, graph);
            std::cout << "before presolving: " << std::endl;
            printResult(grid, active, observed);
            pds::presolve(grid, zero_injection, active, observed);
            std::cout << "after presolving: " << std::endl;
            printResult(grid, active, observed);
        }
        auto active = boost::compose_property_map(boost::associative_property_map(active_data), get(&pds::Bus::id, graph));
        pds::solve_pds(grid, active);
        pds::map<long, bool> observed_map;
        auto observed = boost::compose_property_map(boost::associative_property_map(observed_map), get(&pds::Bus::id, graph));
        pds::dominate(graph, active, observed);
        pds::propagate(graph, get(&pds::Bus::zero_injection, graph),observed);
        printResult(graph, active, observed);
        pds::PdsState state(graph);
        for (auto v: state.vertices()) {
            if (active_data[state[v].id] == pds::PmuState::Active) {
                state.set_active(v);
            }
        }
        for (auto v: state.vertices()) {
            if (!state.is_observed(v)) {
                std::cout << "unobserved " << state.vertex_id(v) << std::endl;
            }
        }
        std::cout << "#unobserved: " << ranges::distance(state.vertices() | ranges::views::filter([&state](auto v) { return state.is_observed(v);})) << std::endl;
        std::cout << "blub " << state.all_observed() << std::endl;
    }
}

int main(int argc, const char** argv) {
    std::string filename = argv[1];
    processBoost(filename);
    return 0;
}
