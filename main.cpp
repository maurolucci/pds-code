#include <iostream>

#include <boost/program_options.hpp>

#include <string>
#include <fmt/format.h>
#include <fmt/chrono.h>
#include <chrono>
#include <functional>

#include <filesystem>

#include "pds.hpp"
#include "graphio.hpp"
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
                    "solver,s",
                    po::value<string>()->default_value("gurobi"),
                    "solve method. Can be any of [none,greedy,greedy-degree,branching,gurobi,brimkov,jovanovic,domination]"
            )
            ("subproblem,u", "split problem into subproblems and solve them individually")
            ("print-solve", "print intermediate solve state")
            ("print-state,p", "print solve state after each step")
            ("time-limit,t", po::value<double>()->default_value(600.0), "time limit for gurobi in seconds")
            ("reductions,r", "apply reductions before exact solving")
            ("export-subproblems,e", "export subproblems as graphml files")
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

    string solve = vm["solver"].as<string>();

    std::function<SolveState(PdsState& state, bool output, double timeLimit)> solver;
    if (solve == "gurobi"s) {
        solver = solve_pds;
    } else if (solve == "none"s) {
        solver = [](auto&, auto, auto) { return SolveState::Other; };
    } else if (solve == "domination"s) {
        solver = solveDominatingSet;
    } else if (solve == "brimkov-orig"s) {
        solver = solveBrimkov;
    } else if (solve == "brimkov"s) {
        solver = solveBrimkovExpanded;
    } else if (solve == "jovanovic"s) {
        solver = solveJovanovic;
    } else if (solve == "greedy"s) {
        solver = [&vm](auto& state, bool, double){ return solveGreedy(state, vm.count("reductions"), greedy_strategies::largestObservationNeighborhood);};
    } else if (solve == "greedy-degree"s) {
        solver = [&vm](auto& state, bool, double){ return solveGreedy(state, vm.count("reductions"), greedy_strategies::largestDegree);};
    } else if (solve == "greedy-median"s) {
        solver = [&vm](auto& state, bool, double){ return solveGreedy(state, vm.count("reductions"), greedy_strategies::medianDegree);};
    } else if (solve == "branching"s) {
        solver = [&vm](auto& state, bool, double){ return solveBranching(state, vm.count("reductions")); };
    } else {
        fmt::print("unrecognized solver {}\n", solve);
    }

    if (!vm.count("graph")) {
        fmt::print(stderr, "no input given\n");
        return 1;
    }

    PdsState state(readAutoGraph(vm["graph"].as<string>(), vm["all-zi"].as<bool>()));
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
        exhaustiveReductions(state, true, drawCallback);
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

    if (vm.count("export-subproblems")) {
        for (size_t i = 0; const auto &sub: subproblems) {
            auto graph = sub.graph();
            auto ids = graph.vertices() | ranges::views::transform([&graph](auto v){ return graph[v].id;}) | ranges::to<set<long>>();
            auto freshId = [&ids]() {
                long id = 0;
                while (ids.contains(id)) {++id;}
                ids.insert(id);
                return id;
            };
            for (auto v: sub.graph().vertices()) {
                if (sub.isActive(v)) {
                    auto v1 = graph.addVertex(Bus{.name="temp", .id=freshId(), .zero_injection=true});
                    auto v2 = graph.addVertex(Bus{.name="temp", .id=freshId(), .zero_injection=true});
                    graph.addEdge(v, v1);
                    graph.addEdge(v, v2);
                }
            }
            std::ofstream outstream(fmt::format("{}/comp_{:03}_unsolved.graphml", outdir, i));
            exportGraphml(graph, outstream);
            ++i;
        }
    }

    SolveState result = SolveState::Optimal;


    if (vm.count("subproblem")) {
        for (size_t i = 0; auto &subproblem: subproblems) {
            auto tSub = now();
            if (vm.count("reductions")) {
                exhaustiveReductions(subproblem, false);
                if (drawOptions.drawSubproblems && drawOptions.drawReductions) {
                    pds::drawGrid(
                            subproblem.graph(), subproblem.active(), subproblem.observed(),
                            fmt::format("{}/comp_{:03}_1reductions.svg", outdir, i), layout);
                }
            }
            if (!subproblem.solveTrivial()) {
                fmt::print("solving subproblem {}\n", i);
                printState(subproblem);
                result = solver(subproblem, vm.count("print-solve"), vm["time-limit"].as<double>());
                for (auto v: subproblem.graph().vertices()) {
                    if (subproblem.isActive(v)) {
                        state.setActive(v);
                    }
                }
                auto tSubEnd = now();
                fmt::print("solved subproblem {} in {}\n", i, ms(tSubEnd - tSub));
                if (drawOptions.drawSubproblems && drawOptions.drawSolution) {
                    pds::drawGrid(
                            subproblem.graph(), subproblem.active(), subproblem.observed(),
                            fmt::format("{}/comp_{:03}_2solved.svg", outdir, i), layout);
                }
                ++i;
            }
        }
    } else {
        result = solver(state, vm.count("print-solve"), vm["time-limit"].as<double>());
    }

    if (result == SolveState::Infeasible) {
        auto tSolveEnd = now();
        fmt::print("model proved infeasible after {}\n", ms(tSolveEnd - tSolveStart));
        return 1;
    } else {
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
        printState(input);

        return 0;
    }
}

int main(int argc, const char** argv) {
    //std::string filename = argv[1];
    //processBoost(filename);
    run(argc, argv);
    return 0;
}
