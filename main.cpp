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
        const pds::PdsState& state
) {
    size_t active_count = 0, inactive_count = 0, zero_injection_count = 0, observed_count = 0;
    auto filter_active = [&state](auto v) -> bool {
        return state.isActive(v);
    };
    if (false) {
        auto active = state.graph().vertices() | ranges::views::filter(filter_active) | ranges::to<std::vector>();
        fmt::print("active nodes: {}\n", active);
    }
    for (const auto &v: state.graph().vertices()) {
        switch (state.activeState(v)) {
            case pds::PmuState::Active:
                ++active_count;
                break;
            case pds::PmuState::Inactive:
                ++inactive_count;
                break;
            case pds::PmuState::Blank:
                break;
        }
        zero_injection_count += state.graph()[v].zero_injection;
        observed_count += state.isObserved(v);
    }
    fmt::print(
                "graph (n={}, m={}, #active={}, #inactive={}, #observed={}, #zero_injection={})\n",
                state.graph().numVertices(),
                state.graph().numEdges(),
                active_count,
                inactive_count,
                observed_count,
                zero_injection_count);
    bool feasible = state.allObserved();
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
            (
                    "all-zi,z",
                    po::value<bool>()->default_value(false)->implicit_value(true, "true"),
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
    } else if (solve == "fast-greedy"s) {
        solver = [](auto& state, bool, double){ return fastGreedy(state);};
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

    if (drawOptions.drawInput) {
        writePds(state.graph(), fmt::format("{}/0_input.pds", outdir));
    }
    auto printState = [&](const PdsState& state) {
        if (vm.count("print-state")) printResult(state);
    };
    fmt::print("input:\n");
    printState(state);

    auto tSolveStart = now();
    size_t counter = 0;

    auto drawCallback = [&](const pds::PdsState &state, const std::string &name) mutable {
        if (drawOptions.drawReductions) {
            writePds(state.graph(), fmt::format("{}/1_red_{:04}_{}.pds", outdir, counter, name) );
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
    ranges::sort(subproblems,
                 [](const pds::PdsState &left, const pds::PdsState &right) -> bool {
                     return left.graph().numVertices() < right.graph().numVertices();
                 });
    if (drawOptions.drawSubproblems) {
        for (size_t i = 0; auto &subproblem: subproblems) {
            if (!subproblem.allObserved()) {
                writePds(subproblem.graph(), fmt::format("{}/comp_{:03}_0unsolved.pds", outdir, i));
                ++i;
            }
        }
    }

    SolveState result = SolveState::Optimal;

    if (vm.count("subproblem")) {
        for (size_t i = 0; auto &subproblem: subproblems) {
            auto tSub = now();
            fmt::print("solving subproblem {}\n", i);
            printState(subproblem);
            result = combineSolveState(result, solver(subproblem, vm.count("print-solve"), vm["time-limit"].as<double>()));
            state.applySubsolution(subproblem);
            auto tSubEnd = now();
            fmt::print("solved subproblem {} in {} ({} active)\n", i, ms(tSubEnd - tSub), subproblem.numActive());
            if (drawOptions.drawSubproblems && drawOptions.drawSolution) {
                writePds(subproblem.graph(), fmt::format("{}/comp_{:03}_2solved.pds", outdir, i));
            }
            ++i;
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
            writePds(state.graph(), fmt::format("{}/2_solved_preprocessed.pds", outdir));
            writePds(input.graph(), fmt::format("{}/3_solved.pds", outdir));
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
