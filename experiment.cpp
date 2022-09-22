#include <chrono>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <fmt/format.h>
#include <fmt/os.h>
#include <fmt/chrono.h>

#include "pds.hpp"
#include "graphio.hpp"
#include "pdssolve.hpp"
#include "gurobi_solve.hpp"

using namespace pds;

// Utility
template<> struct fmt::formatter<SolveState>: formatter<string_view> {
    template<class FormatContext>
    auto format(SolveState state, FormatContext& ctx) {
        string_view name = "unknown";
        switch (state) {
            case pds::SolveState::Optimal: name = "Optimal"; break;
            case pds::SolveState::Other: name = "Other"; break;
            case pds::SolveState::Heuristic: name = "Heuristic"; break;
            case pds::SolveState::Infeasible: name = "Infeasible"; break;
            case pds::SolveState::Timeout: name = "Timeout"; break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

auto now() {
    return std::chrono::high_resolution_clock::now();
}

template<typename T>
auto ms(T time) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(time).count();
}

// Reductions

bool simpleReductions(PdsState& state) {
    return exhaustiveSimpleReductions(state);
}

bool applyDominationReductions(PdsState& state) {
    return dominationReductions(state);
}

bool applyReductions(PdsState& state) {
    return exhaustiveReductions(state);
}

// Main

int main(int argc, const char** argv) {
    namespace po = boost::program_options;
    using namespace std::string_literals;
    using std::string;

    po::options_description desc(argv[0]);
    desc.add_options()
            ("help,h", "show this help")
            ("graph,f", po::value<std::vector<string>>()->required()->multitoken(), "input files")
            ("ignore,i", po::value<string>(), "file that lists inputs to ignore")
            ("all-zi,z", "consider all nodes zero-innjection")
            ("outfile,o", po::value<string>(), "output file")
            ("reduce,r", po::value<string>()->implicit_value("all"s,"all")->default_value("none"s,"none"), "apply reduce. can be any of [none,all,simple,domination]")
            ("solve,s", po::value<string>()->default_value("none")->implicit_value("gurobi"), "solve method. can be any of [none,gurobi,greedy]")
            ("subproblems,u", "split into subproblems before calling solve")
            ("timeout,t", po::value<double>()->default_value(600.), "gurobi time limit (seconds)")
    ;
    po::positional_options_description pos;
    pos.add("graph", -1);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(),vm);

    if (vm.count("help")) {
        desc.print(std::cout);
        return 1;
    }

    bool allZeroInjection = vm.count("all-zi");
    double timeout = vm["timeout"].as<double>();

    string reductionName = vm["reduce"].as<string>();
    std::function<bool(PdsState&)> reduce;
    if (reductionName == "none") {
        reduce = [](auto&) { return false; };
    } else if (reductionName == "all") {
        reduce = applyReductions;
    } else if (reductionName == "simple") {
        reduce = simpleReductions;
    } else if (reductionName == "domination") {
        reduce = applyDominationReductions;
    } else {
        fmt::print(stderr, "invalid reduction mode: {}", reductionName);
        return 2;
    }

    string solverName = vm["solve"].as<string>();
    std::function<SolveState(PdsState&)> solve;
    if (solverName == "gurobi") {
        solve = [timeout](auto &state) { return solve_pds(state, false, timeout); };
    } else if (solverName == "brimkov") {
        solve = [timeout](auto& state) { return solveBrimkovExpanded(state, false, timeout); };
    } else if (solverName == "greedy") {
        solve = [](auto& state) { return solveGreedy(state, true); };
    } else if (solverName == "none") {
        solve = [](auto&) { return SolveState::Other; };
    } else {
        fmt::print(stderr, "invalid solve: {}", solverName);
        return 2;
    }

    bool subproblems = vm.count("subproblems");

    std::vector<string> inputs;

    if (vm.count("graph")) {
        inputs = vm["graph"].as<std::vector<string>>();
    } else {
        fmt::print(stderr, "no input");
        return 2;
    }

    FILE* outfile = nullptr;
    std::vector<FILE*> outputs = {stdout};
    if (vm.count("outfile")) {
        outfile = fopen(vm["outfile"].as<string>().c_str(), "w");
        outputs.push_back(outfile);
    }
    for (auto out: outputs) {
        fmt::print(out, "#{}\n", fmt::join(std::span(argv, argc + argv), " "));
        fmt::print(out, "{}\n","name,pmus,solved,result,t_total,t_reductions,t_solver,n,m,zi,n_reduced,m_reduced,zi_reduced,pmu_reduced,blank_reduced");
    }


    for (const std::string& filename: inputs) {
        PdsState state(readAutoGraph(filename, allZeroInjection));
        size_t n = state.graph().numVertices();
        size_t m = state.graph().numEdges();
        size_t zi = state.numZeroInjection();
        auto t0 = now();
        reduce(state);
        size_t nReduced = state.graph().numVertices();
        size_t mReduced = state.graph().numEdges();
        size_t ziReduced = state.numZeroInjection();
        size_t pmusReduced = state.numActive();
        size_t blankReduced = nReduced - state.numActive() - state.numInactive();
        auto t1 = now();
        SolveState result = SolveState::Optimal;
        if (subproblems) {
            for (auto substate: state.subproblems(true)) {
                result = combineSolveState(result, solve(substate));
            }
        } else {
            result = solve(state);
        }
        auto t2 = now();
        size_t pmus = state.numActive();
        fmt::memory_buffer buf;
        using namespace fmt::literals;
        for (auto file: outputs) {
            fmt::print(file,
                    "{name},{pmus},{solved},{result},{t_total},{t_reductions},{t_solver},{n},{m},{zi},{nReduced},{mReduced},{ziReduced},{pmuReduced},{blankReduced}\n",
                    "name"_a = filename,
                    "pmus"_a = pmus,
                    "solved"_a = state.allObserved(),
                    "result"_a = result,
                    "t_total"_a = ms(t2 - t0),
                    "t_reductions"_a = ms(t1 - t0),
                    "t_solver"_a = ms(t2 - t1),
                    "n"_a = n,
                    "m"_a = m,
                    "zi"_a = zi,
                    "nReduced"_a = nReduced,
                    "mReduced"_a = mReduced,
                    "ziReduced"_a = ziReduced,
                    "pmuReduced"_a = pmusReduced,
                    "blankReduced"_a = blankReduced
            );
        }
    }
    if (outfile != nullptr) {
        fclose(outfile);
    }
}