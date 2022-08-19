#include <iostream>

#include <boost/program_options.hpp>

#include <string>
#include <fmt/format.h>
#include <fmt/chrono.h>
#include <chrono>
#include <functional>

#include <filesystem>

#include "pds.hpp"
#include "pdssolve.hpp"
#include "gurobi_solve.hpp"
#include "draw_grid.hpp"

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

bool simpleReductions(pds::PdsState& state, const std::function<void(const pds::PdsState&, std::string)>& callback) {
    bool changed = false;
    if (state.disableLowDegree()) { callback(state, "low_degree"); changed = true; }
    while (state.collapseLeaves()) { callback(state, "leaves"); changed = true; }
    if (state.reduceObservedNonZi()) { callback(state, "non_zi"); changed = true; }
    while (state.collapseDegreeTwo()) { callback(state, "path"); changed = true; }
    if (state.collapseObservedEdges()) { callback(state, "observed_edges"); changed = true; }
    return changed;
}

bool applyReductions(pds::PdsState& state, const std::function<void(const pds::PdsState&, std::string)>& callback) {
    bool changed = false;
    while (simpleReductions(state, callback)) { changed = true; }
    if (state.disableObservationNeighborhood()) { callback(state, "observation_neighborhood"); changed = true; }
    if (state.activateNecessaryNodes()) { callback(state, "necessary_nodes"); changed = true; }
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

struct DrawOptions {
    bool drawInput;
    bool drawSolution;
    bool drawReductions;
    bool drawSubproblems;

    bool drawAny() { return drawInput || drawSolution || drawReductions || drawSubproblems; }
};

int run(int argc, const char** argv) {
    using namespace pds;
    using std::string, std::vector;
    using namespace std::string_literals;
    namespace po = boost::program_options;
    namespace fs = std::filesystem;
    po::options_description desc("options");
    desc.add_options()
            ("help,h", "show this help")
            ("graph,f", po::value<string>(), "input graph")
            ("outdir,o", po::value<string>()->default_value("out"), "output directory")
            (
                    "solve",
                    po::value<string>()->default_value("subproblem"),
                    "gurobi solve method. Can be any of [none,greedy,full,subproblem]"
            )
            ("print-solve", "print intermediate solve state")
            ("print-state,p", "print solve state after each step")
            ("time-limit,t", po::value<double>()->default_value(600.0), "time limit for gurobi in seconds")
            ("reductions,r", "apply reductions before exact solving")
            (
                    "all-zi,z",
                    po::value<bool>()->default_value(false)->implicit_value(true),
                    "consider all vertices zero-inection"
            )
            (
                    "draw,d",
                    po::value<vector<string>>()->default_value({"none"s}, "none")->implicit_value({"all"s}, "all")->composing(),
                    "can be one of [none,all,input,solution,reductions,subproblems]"
            )
            ;
    po::positional_options_description pos;
    pos.add("graph", 1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        desc.print(std::cout);
        return 1;
    }

    auto outdir = vm["outdir"].as<string>();
    switch (fs::status(outdir).type()) {
        case fs::file_type::none:
        case fs::file_type::not_found:
            fs::create_directories(outdir);
            break;
        case fs::file_type::directory:
        case fs::file_type::symlink:
            break;
        default:
            fmt::print(stderr, "could not create output directory '{}'", outdir);
            return 1;
    }

    DrawOptions drawOptions{};
    auto drawOption = vm["draw"].as<vector<string>>();
    for (auto d: drawOption) {
        if (d == "none") {
            drawOptions.drawSubproblems = drawOptions.drawReductions = drawOptions.drawSolution = drawOptions.drawInput = false;
        } else if (d == "all") {
            drawOptions.drawSubproblems = drawOptions.drawReductions = drawOptions.drawSolution = drawOptions.drawInput = true;
        } else if (d == "input") {
            drawOptions.drawInput = true;
        } else if (d == "solution") {
            drawOptions.drawSolution = true;
        } else if (d == "reductions") {
            drawOptions.drawReductions = true;
        } else if (d == "subproblems") {
            drawOptions.drawSubproblems = true;
        } else {
            fmt::print(stderr, "invalid draw option {}\n", d);
            return 1;
        }
    }

    if (!set<string>{"none"s, "greedy"s, "full"s, "subproblem"s}.contains(vm["solve"].as<string>())) {
        fmt::print(stderr, "invalid solve option: {}\n", vm["solve"].as<string>());
        return 1;
    }

    if (!vm.count("graph")) {
        fmt::print(stderr, "no input given\n");
        return 1;
    }

    PdsState state(import_graphml(vm["graph"].as<string>(), vm["all-zi"].as<bool>()));
    auto input = state;

    map<PowerGrid::vertex_descriptor, Coordinate> layout;
    if (drawOptions.drawAny()) {
        fmt::print("computing layout\n");
        auto t0 = now();
        layout = pds::layoutGraph(state.graph());
        auto t1 = now();
        fmt::print("computed layout in {}\n", ms(t1-t0));
    }

    if (drawOptions.drawInput) {
        drawGrid(state.graph(), state.active(), state.observed(), fmt::format("{}/0_input.svg", outdir), layout);
    }
    auto printState = [&](const PdsState& state) {
        if (vm.count("print-state")) printResult(state.graph(), state.active(), state.observed());
    };
    fmt::print("input:\n");
    printState(state);

    auto tSolveStart = now();
    size_t counter = 0;

    auto drawCallback = [&](const pds::PdsState &state, const std::string &name) mutable {
        if (drawOptions.drawReductions) {
            pds::drawGrid(state.graph(),
                          state.active(),
                          state.observed(),
                          fmt::format("{}/1_red_{:04}_{}.svg", outdir, counter, name),
                          layout);
            ++counter;
        }
    };

    if (vm.count("reductions")) {

        fmt::print("applying reductions\n");
        while (applyReductions(state, drawCallback)) { };
        auto tReductions = now();
        printState(state);
        fmt::print("reductions took {}\n", ms(tReductions - tSolveStart));
    }

    vector subproblems = state.subproblems(true);
    if (drawOptions.drawSubproblems) {
        ranges::sort(subproblems,
                     [](const pds::PdsState &left, const pds::PdsState &right) -> bool {
                         return left.graph().numVertices() < right.graph().numVertices();
                     });
        for (size_t i = 0; auto &subproblem: subproblems) {
            if (!subproblem.allObserved()) {
                pds::drawGrid(
                        subproblem.graph(), subproblem.active(), subproblem.observed(),
                        fmt::format("{}/comp_{:03}_0unsolved.svg", outdir, i), layout);
                ++i;
            }
        }
    }

    for (size_t i = 0; auto &subproblem: subproblems) {
        auto tSub = now();
        if (vm.count("reductions")) {
            applyReductions(subproblem, [](const auto&, const auto&) {});
            if (drawOptions.drawSubproblems && drawOptions.drawReductions) {
                pds::drawGrid(
                        subproblem.graph(), subproblem.active(), subproblem.observed(),
                        fmt::format("{}/comp_{:03}_1reductions.svg", outdir, i), layout);
            }
        }
        if (!subproblem.solveTrivial()) {
            fmt::print("solving subproblem {}\n", i);
            printState(subproblem);
            if (vm["solve"].as<string>() == "subproblem") {
                solve_pds(subproblem, vm.count("print-solve"), vm["time-limit"].as<double>());
                for (auto v: subproblem.graph().vertices()) {
                    if (subproblem.isActive(v)) {
                        state.setActive(v);
                    }
                }
                auto tSubEnd = now();
                fmt::print("solved subproblem {} in {}\n", i, ms(tSubEnd - tSub));
            }
            if (drawOptions.drawSubproblems && drawOptions.drawSolution) {
                pds::drawGrid(
                        subproblem.graph(), subproblem.active(), subproblem.observed(),
                        fmt::format("{}/comp_{:03}_2solved.svg", outdir, i), layout);
            }
            ++i;
        }
    }

    if (vm["solve"].as<string>() == "full") {
        solve_pds(state, vm.count("print-solve"), vm["time-limit"].as<double>());
    } else if (vm["solve"].as<string>() == "greedy") {
        solveGreedy(state);
    }
    for (auto v: state.graph().vertices()) {
        if (state.isActive(v)) {
            input.setActive(v);
        }
    }
    auto tSolveEnd = now();
    if (drawOptions.drawSolution) {
        drawGrid(state.graph(), state.active(), state.observed(), fmt::format("{}/2_solved_preprocessed.svg", outdir), layout);
        drawGrid(input.graph(), input.active(), input.observed(), fmt::format("{}/3_solved.svg", outdir), layout);
    }
    fmt::print("solved in {}\n", ms(tSolveEnd - tSolveStart));
    printState(state);

    return 0;
}

int main(int argc, const char** argv) {
    //std::string filename = argv[1];
    //processBoost(filename);
    run(argc, argv);
    return 0;
}
