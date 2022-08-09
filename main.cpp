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
    size_t active_count = 0, inactive_count = 0, zero_injection_count = 0, observed_count = 0;
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
        observed_count += observed.contains(v);
    }
    std::cout << "]\n";
    std::cout << "graph (m=" << graph.numEdges() << ", n=" << graph.numVertices()
              << ", #active=" << active_count
              << ", #inactive=" << inactive_count
              << ", #observed=" << observed_count
              << ", #zero_injection=" << zero_injection_count << ")\n";
    bool feasible = pds::observed(graph, observed);
    std::cout << "feasible: " << std::boolalpha << feasible << std::endl;
    if (!feasible) {
        bool sep = false;
        std::cout << "active: [";
        for (const auto& v: graph.vertices()) {
            if (active.at(v) == pds::PmuState::Active) {
                if (sep) std::cout << ", ";
                sep = true;
                std::cout << graph[v].name;
            }
        }
        std::cout << "]\n";
        sep = false;
        std::cout << "unobserved: [";
        for (auto v: graph.vertices()) {
            if (!observed.contains(v)) {
                if (sep) std::cout << ", ";
                sep = true;
                std::cout << graph[v].name;
            }
        }
        std::cout << "]\n";
    }
}

void printGraph(
        const pds::PowerGrid& graph
) {
    std::cout << "graph {\n";
    for (auto v: graph.vertices()) {
        std::cout << graph[v].name << "; ";
    }
    std::cout << "\n";
    for (auto e: graph.edges()) {
        std::cout << graph[graph.source(e)].name << " -- " << graph[graph.target(e)].name << ";\n";
    }
    std::cout << "}\n";
}

void processBoost(const std::string& filename) {
    auto graph = pds::import_graphml(filename);
    printGraph(graph);
    for (auto v: graph.vertices()) {
        graph[v].zero_injection = true;//zero_injection[v];
    }
    {
        auto active = graph.vertices()
        | ranges::views::transform([](pds::PowerGrid::vertex_descriptor v) { return std::make_pair(v, pds::PmuState::Blank);})
        | ranges::to<pds::map<pds::PowerGrid::vertex_descriptor, pds::PmuState>>();
        pds::solve_pds(graph, active);
        pds::set<pds::PowerGrid::vertex_descriptor> observed;
        printResult(graph, active, observed);
        pds::dominate(graph, active, observed);
        pds::propagate(graph, observed);
        printResult(graph, active, observed);
        pds::PdsState state(graph);
        for (auto v: graph.vertices()) {
            if (active.at(v) == pds::PmuState::Active) {
                state.setActive(v);
            }
        }
        std::cout << "#unobserved: " << ranges::distance(state.graph().vertices() | ranges::views::filter([&state](auto v) { return !state.is_observed(v);})) << std::endl;
    }
    {
        pds::PdsState state(graph);
        while (state.collapseLeaves()) { }
        while (state.collapseDegreeTwo()) { }
        auto active = graph.vertices()
                      | ranges::views::transform([](pds::PowerGrid::vertex_descriptor v) { return std::make_pair(v, pds::PmuState::Blank);})
                      | ranges::to<pds::map<pds::PowerGrid::vertex_descriptor, pds::PmuState>>();
        //try {
            pds::solve_pds(state.graph(), active);
        //} catch (GRBException ex) {
        //    std::cerr << "gurobi exception: " << ex.getErrorCode() << " (" << ex.getMessage() << ")\n";
        //    throw ex;
        //}
        pds::set<pds::PowerGrid::vertex_descriptor> observed;
        pds::dominate(graph, active, observed);
        pds::propagate(graph, observed);
        printResult(graph, active, observed);
        std::cout << "#unobserved: " << ranges::distance(state.graph().vertices() | ranges::views::filter([&state](auto v) { return !state.is_observed(v);})) << std::endl;
    }
}

int main(int argc, const char** argv) {
    std::string filename = argv[1];
    processBoost(filename);
    return 0;
}
