#include <iostream>

#include <boost/graph/graphml.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/property_map/compose_property_map.hpp>

#include <string>
#include <fmt/format.h>
#include <fmt/chrono.h>
#include <chrono>
#include <functional>

#include "pds.hpp"
#include "gurobi_solve.hpp"
#include "draw_grid.hpp"

template<class T, class Graph>
auto make_vector_property_map(T&& vec, Graph graph) {
    return boost::make_iterator_property_map(vec.begin(), get(boost::vertex_index, graph));
}

void printResult(
        const pds::PowerGrid& graph,
        const pds::map<pds::PowerGrid::vertex_descriptor, pds::PmuState>& active,
        const pds::set<pds::PowerGrid::vertex_descriptor>& observed
) {
    size_t active_count = 0, inactive_count = 0, zero_injection_count = 0, observed_count = 0;
    auto filter_active = [&active](auto v) -> bool {
        return active.at(v) == pds::PmuState::Active;
    };
    if (false) {
        auto active = graph.vertices() | ranges::views::filter(filter_active) | ranges::to<std::vector>();
        fmt::print("active nodes: {}\n", active);
    }
    for (const auto &v: graph.vertices()) {
        switch (active.at(v)) {
            case pds::PmuState::Active:
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
    fmt::print(
                "graph (n={}, m={}, #active={}, #inactive={}, #observed={}, #zero_injection={})\n",
                graph.numVertices(),
                graph.numEdges(),
                active_count,
                inactive_count,
                observed_count,
                zero_injection_count);
    bool feasible = pds::observed(graph, observed);
    fmt::print("solved: {}\n", feasible);
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

auto now() {
    return std::chrono::high_resolution_clock::now();
}

template<typename T>
auto ms(T time) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(time);
}

bool applyReductions(pds::PdsState& state, const std::function<void(const pds::PdsState&, std::string)>& callback) {
    bool changed = false;
    while (state.collapseLeaves()) { callback(state, "leaves"); changed = true; }
    if (state.reduceObservedNonZi()) { callback(state, "non_zi"); changed = true; }
    while (state.disableObservationNeighborhood()) { callback(state, "observation_neighborhood"); changed = true; }
    while (state.activateNecessaryNodes()) { callback(state, "necessary_nodes"); changed = true; }
    while (state.collapseDegreeTwo()) { callback(state, "path"); changed = true; }
    if (state.collapseObservedEdges()) { callback(state, "observed_edges"); changed = true; }
    return changed;
}

void processBoost(const std::string& filename) {
    fmt::print("loading graph...\n");
    auto t0 = now();
    auto graph = pds::import_graphml(filename, true);
    auto t1 = now();
    fmt::print("graph loaded in {}\n", ms(t1 - t0)); std::fflush(nullptr);
    //printGraph(graph);
    fmt::print("computing layout...\n");
    auto layout = pds::layoutGraph(graph);
    auto t2 = now();
    fmt::print("layout computed in {}\n", ms(t2 - t1));
    if (graph.numVertices() < 1000) {
        fmt::print("solve exact\n");
        auto t3 = now();
        auto active = graph.vertices()
        | ranges::views::transform([](pds::PowerGrid::vertex_descriptor v) { return std::make_pair(v, pds::PmuState::Blank);})
        | ranges::to<pds::map<pds::PowerGrid::vertex_descriptor, pds::PmuState>>();
        pds::solve_pds(graph, active);
        auto t4 = now();
        fmt::print("solved in {}", ms(t4 - t3));
        pds::set<pds::PowerGrid::vertex_descriptor> observed;
        printResult(graph, active, observed);
        pds::dominate(graph, active, observed);
        pds::propagate(graph, observed);
        printResult(graph, active, observed);
        pds::drawGrid(graph, active, observed, "out/4_solve_input.svg", layout);
    }
    {
        fmt::print("solve with reductions\n");
        pds::PdsState state(graph);
        bool drawReductionSteps = true;
        size_t counter = 0;
        auto drawState = [&](const pds::PdsState& state, const std::string& name) mutable {
            if (drawReductionSteps) {
                pds::drawGrid(state.graph(), state.active(), state.observed(), fmt::format("out/1_red_{:04}_{}.svg", counter, name), layout);
                ++counter;
            }
        };
        pds::drawGrid(state.graph(), state.active(), state.observed(), "out/0_input.svg", layout);
        fmt::print("applying reductions\n");
        auto t_startReduction = now();
        while (applyReductions(state, drawState)) { }
        auto t_endReduction = now();
        fmt::print("reductions took {}\n", ms(t_endReduction - t_startReduction));
        printResult(state.graph(), state.active(), state.observed());
        pds::drawGrid(state.graph(), state.active(), state.observed(), "out/1_red_preprocessed.svg", layout);
        bool drawTrivial = false;
        auto t_solveStart = now();
        bool optimal = true;
        if (true) {
            auto subproblems = state.subproblems(true);
            ranges::sort(subproblems,
                         [](const pds::PdsState &left, const pds::PdsState &right) -> bool {
                             return left.graph().numVertices() < right.graph().numVertices();
                         });
            for (size_t i = 0, j = 0; auto &subproblem: subproblems) {
                if (subproblem.solveTrivial()) {
                    if (drawTrivial) {
                        pds::drawGrid(subproblem.graph(), subproblem.active(), subproblem.observed(),
                                      fmt::format("out/comp_trivial_{:03}.svg", i), layout);
                        ++i;
                    }
                } else {
                    fmt::print("solving subproblem {}\n", j);
                    printResult(subproblem.graph(), subproblem.active(), subproblem.observed());
                    pds::drawGrid(subproblem.graph(), subproblem.active(), subproblem.observed(),
                                  fmt::format("out/comp_{:03}_0unsolved.svg", j), layout);
                    while (applyReductions(subproblem, [](const auto&, const auto&) { }));
                    pds::drawGrid(subproblem.graph(), subproblem.active(), subproblem.observed(),
                                  fmt::format("out/comp_{:03}_1reductions.svg", j), layout);
                    auto active = subproblem.active();
                    auto t_subSolve = now();
                    if (!pds::solve_pds(subproblem.graph(), active, false, 60 * 60)) {
                        fmt::print("infeasible subproblem {}\n", j);
                    }
                    auto t_subEnd = now();
                    for (auto v: subproblem.graph().vertices()) {
                        if (active.at(v) == pds::PmuState::Active) {
                            subproblem.setActive(v);
                        }
                    }
                    fmt::print("subproblem solved in {}\n", ms(t_subEnd - t_subSolve));
                    fmt::print("solve result:\n");
                    printResult(subproblem.graph(), subproblem.active(), subproblem.observed());
                    pds::drawGrid(subproblem.graph(), subproblem.active(), subproblem.observed(),
                                  fmt::format("out/comp_{:03}_2solved.svg", j), layout);
                    ++j;
                }
                for (auto v: subproblem.graph().vertices()) {
                    if (subproblem.isActive(v)) {
                        state.setActive(v);
                    }
                }
            }
        } else {
            auto active = state.active();
            if (!pds::solve_pds(state.graph(), active)) {
                fmt::print("infeasible\n");
            }
            for (auto v: state.graph().vertices()) {
                if (active.at(v) == pds::PmuState::Active) {
                    state.setActive(v);
                }
            }
        }
        auto t_solveEnd = now();
        if (state.allObserved()) {
            fmt::print("solved in {}\n", ms(t_solveEnd - t_solveStart));
        } else {
            fmt::print("not solved in {}\n", ms(t_solveEnd - t_solveStart));
        }
        pds::drawGrid(state.graph(), state.active(), state.observed(), "out/2_solved_preprocessed.svg", layout);
        auto active = state.active();
        pds::set<pds::PowerGrid::vertex_descriptor> observed;
        pds::dominate(graph, active, observed);
        pds::propagate(graph, observed);
        pds::drawGrid(graph, active, observed, "out/3_solved.svg", layout);
        printResult(graph, active, observed);
        auto isUnobserved = [&state] (auto v) -> bool { return !state.isObserved(v); };
        fmt::print("#unobserved = {}", ranges::distance(state.graph().vertices() | ranges::views::filter(isUnobserved)));
    }
}

int main(int argc, const char** argv) {
    std::string filename = argv[1];
    processBoost(filename);
    return 0;
}
