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

void printResult(
        const pds::PowerGrid& graph,
        const pds::map<pds::PowerGrid::vertex_descriptor, pds::PmuState>& active,
        const pds::set<pds::PowerGrid::vertex_descriptor>& observed
) {
    bool sep = false;
    size_t active_count = 0, inactive_count = 0, zero_injection_count = 0;
    std::cout << "active nodes: [";
    for (const auto &v: graph.vertices()) {
        switch (active.at(v)) {
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
    std::cout << "graph (m=" << graph.numEdges() << ", n=" << graph.numVertices()
              << ", #active=" << active_count
              << ", #inactive=" << inactive_count
              << ", #zero_injection=" << zero_injection_count << ")\n";
    std::cout << "feasible: " << std::boolalpha << pds::observed(graph, observed) << std::endl;
}

void processBoost(const std::string& filename) {
    auto graph = pds::import_graphml(filename);
    for (auto v: graph.vertices()) {
        graph[v].zero_injection = true;//zero_injection[v];
    }
    {
        pds::PowerGrid grid(graph);
        auto active = graph.vertices()
        | ranges::views::transform([](pds::PowerGrid::vertex_descriptor v) { return std::make_pair(v, pds::PmuState::Blank);})
        | ranges::to<pds::map<pds::PowerGrid::vertex_descriptor, pds::PmuState>>();
        //{
        //    pds::set<pds::PowerGrid::vertex_descriptor> observed_data;
        //    std::cout << "before presolving: " << std::endl;
        //    printResult(grid, active_data, observed_data);
        //    pds::presolve(grid, zero_injection, active, observed);
        //    std::cout << "after presolving: " << std::endl;
        //    printResult(grid, active, observed);
        //}
        pds::solve_pds(grid, active);
        pds::set<pds::PowerGrid::vertex_descriptor> observed;
        printResult(graph, active, observed);
        pds::dominate(graph, active, observed);
        pds::propagate(graph,observed);
        printResult(graph, active, observed);
        pds::PdsState state(graph);
        std::cout << "#unobserved: " << ranges::distance(state.graph().vertices() | ranges::views::filter([&state](auto v) { return !state.is_observed(v);})) << std::endl;
        std::cout << "blub " << state.all_observed() << std::endl;
    }
}

int main(int argc, const char** argv) {
    std::string filename = argv[1];
    processBoost(filename);
    return 0;
}
