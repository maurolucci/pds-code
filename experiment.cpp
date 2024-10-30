#include <chrono>
#include <string>
#include <vector>
#include <optional>

#include <boost/program_options.hpp>
#include <fmt/format.h>
#include <fmt/chrono.h>

#include <htd/main.hpp>

#include "pds.hpp"
#include "graphio.hpp"
#include "pdssolve.hpp"
#include "gurobi_solve.hpp"
#include "fort_solve.hpp"
#include "cycle_solve.hpp"
#include "cycle_fort_solve.hpp"
#include "path_solve.hpp"
#include "path_fort_solve.hpp"

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
auto µs(T time) {
    return std::chrono::duration_cast<std::chrono::microseconds>(time).count();
}

// Reductions

bool simpleReductions(PdsState& state) {
    return exhaustiveSimpleReductions(state);
}

bool applyDominationReductions(PdsState& state) {
    bool changed = false;
    while (dominationReductions(state)) { changed = true; }
    return changed;
}

bool applyReductionsNotDomination(PdsState& state) {
    bool anyChanged = false;
    bool firstRun = true;
    bool changed;
    do {
        changed = exhaustiveSimpleReductions(state);
        if (firstRun || changed) if(state.activateNecessaryNodes()) {
                changed = true;
            }
        firstRun = false;
        anyChanged |= changed;
    } while (changed);
    return anyChanged;
}

bool applyReductionsNotNecessary(PdsState& state) {
    return noNecessaryReductions(state);
}

bool applyReductions(PdsState& state) {
    return exhaustiveReductions(state);
}

size_t treeWidth(const PowerGrid& graph) {
    std::unique_ptr<htd::LibraryInstance> library(htd::createManagementInstance(htd::Id::FIRST));
    htd::Graph htdGraph(library.get());
    pds::map<PowerGrid::VertexDescriptor, htd::vertex_t> vertices;
    for (auto v: graph.vertices()){
        auto mapped = htdGraph.addVertex();
        vertices[v] = mapped;
    }
    for (auto edge: graph.edges()) {
        auto [s, t] = graph.endpoints(edge);
        htdGraph.addEdge(vertices[s], vertices[t]);
    }
    library->orderingAlgorithmFactory().setConstructionTemplate(new htd::MinFillOrderingAlgorithm(library.get()));
    auto treeDecompositionAlgo = std::make_unique<htd::CombinedWidthMinimizingTreeDecompositionAlgorithm>(library.get());
    treeDecompositionAlgo->setComputeInducedEdgesEnabled(false);
    auto algo = std::make_unique<htd::CombinedWidthMinimizingTreeDecompositionAlgorithm>(library.get());
    auto baseAlgo = std::make_unique<htd::WidthMinimizingTreeDecompositionAlgorithm>(library.get());
    baseAlgo->setIterationCount(10);
    baseAlgo->addManipulationOperation(new htd::NormalizationOperation(library.get()));
    algo->addDecompositionAlgorithm(baseAlgo.release());
    auto decomposition = std::unique_ptr<htd::ITreeDecomposition>(algo->computeDecomposition(htdGraph));
    return decomposition->maximumBagSize();
}

VertexMap<size_t> propagationDistance(const PowerGrid& graph)  {
    using Vertex = PowerGrid::VertexDescriptor;
    VertexMap<ssize_t> unobservedDegree;
    VertexMap<size_t> step;
    std::deque<Vertex> queue;
    for (auto v: graph.vertices()) {
        unobservedDegree[v] += graph.degree(v);
        if (graph.getVertex(v).pmu == PmuState::Active) {
            step.insert_or_assign(v, 0);
            for (auto w: graph.neighbors(v)) {
                unobservedDegree[w] -= 1;
                if (!step.contains(w) && graph.getVertex(w).pmu != pds::PmuState::Active) {
                    step.insert_or_assign(w, step.at(v) + 1);
                    for (auto u: graph.neighbors(w)) {
                        unobservedDegree[u] -= 1;
                    }
                }
            }
        }
    }
    for (auto v: graph.vertices()) {
        if (step.contains(v) && unobservedDegree[v] == 1) queue.push_back(v);
    }
    for (auto v: graph.vertices()) {
        auto neighbors = graph.neighbors(v) | ranges::to<std::vector>;
        size_t unobserved = std::ranges::distance(neighbors | ranges::views::filter([&step](auto v) { return !step.contains(v); }));
        unused(unobserved);
        assert(unobservedDegree[v] == unobserved);
        assert(!step.contains(v) || step[v] == (graph.getVertex(v).pmu != PmuState::Active));
        if (graph.getVertex(v).pmu == pds::PmuState::Active) {
            for (auto w: graph.neighbors(v)) {
                assert(step.contains(w));
            }
        }
    }
    while (!queue.empty()) {
        auto v = queue.front();
        queue.pop_front();
        if (unobservedDegree[v] == 1) {
            for (auto w: graph.neighbors(v)) {
                if (!step.contains(w)) {
                    step[w] = step[v] + 1;
                    for (auto u: graph.neighbors(w)) {
                        unobservedDegree[u] -= 1;
                        if (unobservedDegree[u] == 1 && step.contains(u)) {
                            queue.push_back(u);
                        }
                    }
                    if (unobservedDegree[w] == 1) {
                        queue.push_back(w);
                    }
                }
            }
        }
    }
    for (auto v: graph.vertices()) {
        assert(unobservedDegree[v] == std::ranges::distance(graph.neighbors(v) | ranges::views::filter([&step](auto v) { return !step.contains(v); })));
    }
    return step;
}

void writeSolutionStatistics(const std::string_view& name, const PdsState& state, FILE* out) {
    using namespace fmt::literals;
    size_t maxDegree = 0;
    map<std::pair<pds::PmuState, bool>, map<size_t, size_t>> degrees;
    for (auto v: state.graph().vertices()) {
        auto deg = state.graph().degree(v);
        maxDegree = std::max(maxDegree, deg);
        degrees[{state.activeState(v), state.isZeroInjection(v)}][deg] += 1;
    }
    std::vector<std::string> degreeCount;
    for (size_t i = 0; i <= maxDegree; ++i) {
        std::vector<size_t> deg;
        for (bool zi: {true, false}) {
            for (auto state: {PmuState::Inactive, PmuState::Blank, PmuState::Active}) {
                deg.push_back(degrees[{state, zi}][i]);
            }
        }
        degreeCount.push_back(fmt::format("{}", fmt::join(deg, ":")));
    }
    auto step = propagationDistance(state.graph());
    ssize_t maxStep = -1;
    for (auto v: step) maxStep = std::max(maxStep, ssize_t(v.second));
    fmt::print(
            out,
            "{name},{n},{m},{n_zero_injection},{n_pmu},{n_inactive},{n_blank},{n_observed},{propagation_distance},{tree_width},\"{degrees}\"\n",
            "name"_a=name,
            "n"_a=state.graph().numVertices(),
            "m"_a=state.graph().numEdges(),
            "n_zero_injection"_a=state.numZeroInjection(),
            "n_pmu"_a=state.numActive(),
            "n_inactive"_a=state.numInactive(),
            "n_blank"_a=state.graph().numVertices() - state.numActive() - state.numInactive(),
            "n_observed"_a=state.numObserved(),
            "propagation_distance"_a=maxStep,
            "tree_width"_a=treeWidth(state.graph()),
            "degrees"_a=fmt::join(degreeCount, ";")
    );
}

struct FortStats {
    double avgSize;
    size_t fortCount;
    size_t numHs;
};

void writeFortStats(const std::string_view& name, size_t run, const FortStats stat, FILE* out) {
    using namespace fmt::literals;
    // name,run,forts,avg_size,num_hitting_sets
    fmt::print(out, "{},{},{},{},{}\n", name, run, stat.fortCount, stat.avgSize, stat.numHs);
}

// Main

auto getModel(const std::string& name) {
    if (name == "gurobi" || name == "jovanovic2") {
        return pds::modelJovanovicExpanded;
    } else if (name == "jovanovic") {
        return pds::modelJovanovic;
    } else if (name == "brimkov") {
        return pds::modelBrimkov;
    } else if (name == "brimkov-geq") {
        return pds::modelBrimkovGeq;
    } else if (name == "brimkov-alt") {
        return pds::modelBrimkovAlt;
    } else if (name == "brimkov2") {
        return pds::modelBrimkovExpanded;
    } else if (name == "azami" || name == "azami-brimkov") {
        return pds::modelAzamiBrimkov;
    } else if (name == "domination") {
        return pds::modelDomination;
    } else {
        throw std::invalid_argument("unknown model " + name);
    }
}
using Solver = std::function<SolveResult(PdsState&, double)>;
namespace po = boost::program_options;
namespace fs = std::filesystem;

auto getSolver(po::variables_map& vm, FortStats& fortStats, FortStats& cycleStat, FortStats& pathStat, const std::string& currentName, size_t run, size_t subproblemNumber, BoundCallback boundCB) {
    using std::string;
    string solverName = vm["solve"].as<string>();
    int verbose = vm.count("verbose") ? vm["verbose"].as<int>() : 0;
    int earlyStop = vm["early-stop"].as<int>();
    int greedyBoundMode = vm["greedy-bounds"].as<int>();
    int fortInit = vm["fort-init"].as<int>();
    int cycleInit = vm["cycle-init"].as<int>();
    int pathInit = vm["path-init"].as<int>();
    int cycleLimit = vm["cycle-limit"].as<int>();
    int pathLimit = vm["path-limit"].as<int>();
    bool intermediateForts = vm.count("intermediate");
    bool intermediateCycles = vm.count("intermediate-cycles");
    bool intermediatePaths = vm.count("intermediate-paths");
    callback::FortCallback fortCallback = [&,solved=size_t{0}]  (callback::When when, const PdsState& state, const std::vector<VertexList>& forts, size_t lower, size_t upper) mutable {
        if (vm.count("fort-stats") && (when == pds::callback::When::FINAL)) {
            size_t totalFortSize = 0;
            for (auto& fort: forts) {
                totalFortSize += fort.size();
            }
            double averageSize = double(totalFortSize) / double(forts.size());
            fortStats = {averageSize, forts.size(), solved};
        }
        ++solved;
    };
    cyclecallback::CycleCallback cycleCallback = [&,solved=size_t{0}]  (cyclecallback::When when, const PdsState& state, const std::vector<VertexList>& cycles, size_t lower, size_t upper) mutable {
        if (vm.count("fort-stats") && when == pds::cyclecallback::When::FINAL) {
            size_t totalCycleSize = 0;
            for (auto& cycle: cycles) {
                totalCycleSize += cycle.size();
            }
            double averageSize = double(totalCycleSize) / double(cycles.size());
            cycleStat = {averageSize, cycles.size(), solved};
        }
        ++solved;
    };
    pathcallback::PathCallback pathCallback = [&,solved=size_t{0}]  (pathcallback::When when, const PdsState& state, const std::vector<VertexList>& paths, size_t lower, size_t upper) mutable {
        if (vm.count("fort-stats") && when == pds::pathcallback::When::FINAL) {
            size_t totalPathSize = 0;
            for (auto& path: paths) {
                totalPathSize += path.size();
            }
            double averageSize = double(totalPathSize) / double(paths.size());
            pathStat = {averageSize, paths.size(), solved};
        }
        ++solved;
    };
    if (solverName == "branching") {
        return Solver{[](auto &state, double) { return solveBranching(state, true, greedy_strategies::largestDegree); }};
    } else if (solverName == "greedy") {
        return  Solver{[](auto &state, double) { return solveGreedy(state, true, greedy_strategies::largestDegree); }};
    } else if (solverName == "fast") {
        return Solver{[](auto &state, double) { return fastGreedy(state, true); }};
    } else if (solverName == "topdown") {
        return Solver{[](auto &state, double) { return topDownGreedy(state); }};
    } else if (solverName == "none") {
        return Solver{[](auto &state, double) {
            return SolveResult{size_t{0}, state.numActive() + state.numBlank(), SolveState::Other};
        }};
    } else if (solverName == "smith") {
        return Solver{[=](auto &state, double timeLimit) {
            state.disableLowDegree();
            return solveBozeman(state, verbose, timeLimit, 1, fortInit,
                                greedyBoundMode, earlyStop, fortCallback, boundCB, -1);
        }};
    } else if (solverName == "bozeman") {
        return Solver{[=](auto &state, double timeLimit) {
            return solveBozeman(state, verbose, timeLimit, 0, fortInit,
                                greedyBoundMode, earlyStop, fortCallback, boundCB, intermediateForts ? 0 : -1);
        }};
    } else if (solverName == "bozeman2") {
        return Solver{[=](auto &state, double timeLimit) {
            state.disableLowDegree();
            return solveBozeman(state, verbose, timeLimit, 2, fortInit,
                                greedyBoundMode, earlyStop, fortCallback, boundCB, intermediateForts ? 2 : -1);
        }};
    } else if (solverName == "bozeman3") {
        return Solver{[=](auto &state, double timeLimit) {
            return solveBozeman(state, verbose, timeLimit, 3, fortInit,
                                greedyBoundMode, earlyStop, fortCallback, boundCB, intermediateForts ? 3 : -1);
        }};
    } else if (solverName == "forts") {
        return Solver{[=](auto &state, double timeLimit) {
            return solveBozeman(state, verbose, timeLimit, 4, fortInit,
                                greedyBoundMode, earlyStop, fortCallback, boundCB, intermediateForts ? 4 : -1);
        }};
    } else if (solverName == "lazy-forts") {
        return Solver{[=](auto &state, double timeLimit) {
            return solveLazyForts(state, verbose, timeLimit, fortCallback, boundCB);
        }};
    } else if (solverName == "cycles") {
        return Solver{[=](auto &state, double timeLimit) {
            return solveCycles(state, verbose, timeLimit, 0, cycleInit,
                               greedyBoundMode, earlyStop, cycleCallback, boundCB, intermediateCycles ? 0 : -1, cycleLimit);
        }};
    } else if (solverName == "cycles-forts") {
        return Solver{[=](auto &state, double timeLimit) {
            return solveCyclesForts(state, verbose, timeLimit, 0, cycleInit,
                                    greedyBoundMode, earlyStop, fortCallback, cycleCallback, boundCB, -1, cycleLimit);
        }};
    } else if (solverName == "cycles2") {
        return Solver{[=](auto &state, double timeLimit) {
            return solveCycles(state, verbose, timeLimit, 1, cycleInit,
                               greedyBoundMode, earlyStop, cycleCallback, boundCB, intermediateCycles ? 0 : -1, cycleLimit);
        }};
    } else if (solverName == "cycles-forts2") {
        return Solver{[=](auto &state, double timeLimit) {
            return solveCyclesForts(state, verbose, timeLimit, 1, cycleInit,
                                    greedyBoundMode, earlyStop, fortCallback, cycleCallback, boundCB, -1, cycleLimit);
        }};
    } else if (solverName == "paths") {
    return Solver{[=](auto &state, double timeLimit) {
        return solvePaths(state, verbose, timeLimit, 0, pathInit,
                            greedyBoundMode, earlyStop, pathCallback, boundCB, -1, pathLimit);
        }};
    } else if (solverName == "paths11") {
    return Solver{[=](auto &state, double timeLimit) {
        return solvePaths(state, verbose, timeLimit, 11, pathInit,
                            greedyBoundMode, earlyStop, pathCallback, boundCB, -1, pathLimit);
        }};
    } else if (solverName == "paths12") {
    return Solver{[=](auto &state, double timeLimit) {
        return solvePaths(state, verbose, timeLimit, 12, pathInit,
                            greedyBoundMode, earlyStop, pathCallback, boundCB, -1, pathLimit);
        }};
    } else if (solverName == "paths13") {
    return Solver{[=](auto &state, double timeLimit) {
        return solvePaths(state, verbose, timeLimit, 13, pathInit,
                            greedyBoundMode, earlyStop, pathCallback, boundCB, -1, pathLimit);
        }};
    } else if (solverName == "paths21") {
    return Solver{[=](auto &state, double timeLimit) {
        return solvePaths(state, verbose, timeLimit, 21, pathInit,
                            greedyBoundMode, earlyStop, pathCallback, boundCB, -1, pathLimit);
        }};
    } else if (solverName == "paths22") {
    return Solver{[=](auto &state, double timeLimit) {
        return solvePaths(state, verbose, timeLimit, 22, pathInit,
                            greedyBoundMode, earlyStop, pathCallback, boundCB, -1, pathLimit);
        }};
    } else if (solverName == "paths31") {
    return Solver{[=](auto &state, double timeLimit) {
        return solvePaths(state, verbose, timeLimit, 31, pathInit,
                            greedyBoundMode, earlyStop, pathCallback, boundCB, -1, pathLimit);
        }};
    } else if (solverName == "paths32") {
    return Solver{[=](auto &state, double timeLimit) {
        return solvePaths(state, verbose, timeLimit, 32, pathInit,
                            greedyBoundMode, earlyStop, pathCallback, boundCB, -1, pathLimit);
        }};
    } else if (solverName == "paths33") {
    return Solver{[=](auto &state, double timeLimit) {
        return solvePaths(state, verbose, timeLimit, 33, pathInit,
                            greedyBoundMode, earlyStop, pathCallback, boundCB, -1, pathLimit);
        }};
    } else if (solverName == "paths34") {
    return Solver{[=](auto &state, double timeLimit) {
        return solvePaths(state, verbose, timeLimit, 34, pathInit,
                            greedyBoundMode, earlyStop, pathCallback, boundCB, -1, pathLimit);
        }};
    } else if (solverName == "paths1321") {
    return Solver{[=](auto &state, double timeLimit) {
        return solvePaths(state, verbose, timeLimit, 1321, pathInit,
                            greedyBoundMode, earlyStop, pathCallback, boundCB, -1, pathLimit);
        }};
    } else if (solverName == "paths-forts") {
        return Solver{[=](auto &state, double timeLimit) {
            return solvePathsForts(state, verbose, timeLimit, 0, pathInit,
                                    greedyBoundMode, earlyStop, fortCallback, pathCallback, boundCB, -1, pathLimit);
        }};     
    } else if (solverName == "paths-forts2") {
        return Solver{[=](auto &state, double timeLimit) {
            return solvePathsForts(state, verbose, timeLimit, 1, pathInit,
                                    greedyBoundMode, earlyStop, fortCallback, pathCallback, boundCB, -1, pathLimit);
        }};
    } else {
        try {
            return Solver{[model = getModel(solverName), verbose,boundCB,fortInit](auto &state, double timeout) {
                auto mip = model(state);
                if (fortInit) {
                    addFortConstraints(mip, state, fortInit);
                }
                auto result = solveMIP(state, mip, verbose, timeout, boundCB, true);
                if (result.state != SolveState::Infeasible) {
                    applySolution(state, mip);
                }
                return result;
            }};
        } catch (std::invalid_argument &ex) {
            fmt::print(stderr, "{}", ex.what());
            throw ex;
        }
    }
}

int main(int argc, const char** argv) {
    namespace po = boost::program_options;
    namespace fs = std::filesystem;
    using std::string;
    using namespace std::string_literals;

    po::options_description desc(argv[0]);
    desc.add_options()
            ("help,h", "show this help")
            ("graph,f", po::value<std::vector<string>>()->required()->multitoken(), "input files")
            ("all-zi,z", "consider all nodes zero-innjection")
            ("outfile,o", po::value<string>(), "output file")
            ("reduce,r", po::value<string>()->implicit_value("all"s,"all")->default_value("none"s,"none"), "apply reduce. can be any of [none,all,simple,domination]")
            ("repeat,n", po::value<size_t>()->default_value(1)->implicit_value(5), "number of experiment repetitions")
            ("solve,s", po::value<string>()->default_value("none")->implicit_value("gurobi"), "solve method. can be any of [none,gurobi,greedy]")
            ("subproblems,u", "split into subproblems before calling solve")
            ("timeout,t", po::value<double>()->default_value(600.), "gurobi time limit (seconds)")
            ("draw,d", po::value<string>()->implicit_value("out"s), "draw states")
            ("write,w", po::value<string>()->implicit_value("solutions"s), "write solutions to the specified directory")
            ("stat-file", po::value<string>(), "write statistics about solutions")
            ("greedy-bounds,b", po::value<int>()->default_value(0)->implicit_value(1), "if possible, use a greedy algorithm to compute an upper bound (0: never, 1: without reductions (faster), 2: with reductions (more precise))")
            ("fort-stats", po::value<string>(), "file for fort statistics")
            ("fort-init,i", po::value<int>()->implicit_value(3)->default_value(0), "fort initialization method")
            ("cycle-init", po::value<int>()->implicit_value(0)->default_value(0), "cycle initialization method")
            ("path-init", po::value<int>()->implicit_value(0)->default_value(0), "path initialization method")
            ("cycle-limit", po::value<int>()->implicit_value(1000000)->default_value(1000000), "maximum number of cycles added in each reoptimization")
            ("path-limit", po::value<int>()->implicit_value(1000000)->default_value(1000000), "maximum number of paths added in each reoptimization")
            ("intermediate", "find violated forts in every new feasbile solution")
            ("intermediate-cycles", "find violated cycles in every new feasbile solution")
            ("intermediate-paths", "find violated paths in every new feasbile solution")
            ("early-stop", po::value<int>()->default_value(0)->implicit_value(2), "stop hitting set solver when violating hitting set is found [0: only optimal, 1: keep lower bound, 2: stop early]")
            ("write-forts", po::value<string>()->implicit_value("hs"), "directory to which to write hitting set instance")
            ("write-bounds", po::value<string>(), "write bound time series")
            ("verbose,v", po::value<int>()->implicit_value(1), "print additional solver status info")
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

    size_t repetitions = vm["repeat"].as<size_t>();

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
    } else if (reductionName == "no-necessary") {
        reduce = applyReductionsNotNecessary;
    } else if (reductionName == "no-domination") {
        reduce = applyReductionsNotDomination;
    } else {
        fmt::print(stderr, "invalid reduction mode: {}. modes: {}", reductionName, fmt::join({"all", "simple", "domination", "no-necessary", "no-domination"}, ", "));
        return 2;
    }

    preloadMIPSolver();
    if (vm.count("write-forts")) {
        std::string fortsDirName = vm["write-forts"].as<string>();
        fs::create_directories(fs::absolute(fs::path(fortsDirName)));
    }
    FortStats fortStats;
    FortStats cycleStats;
    FortStats pathStats;

    std::optional<fs::path> drawdir;

    if (vm.count("draw")) {
        drawdir = vm["draw"].as<string>();
        if (!fs::is_directory(*drawdir)) {
            fs::create_directories(*drawdir);
        }
    }
    std::optional<fs::path> solDir;
    if (vm.count("write")) {
        solDir = vm["write"].as<string>();
        if (!fs::is_directory(*solDir)) {
            fs::create_directories(*solDir);
        }
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
        auto outfileName = vm["outfile"].as<string>();
        fs::create_directories(fs::absolute(fs::path(outfileName)).parent_path());
        outfile = fopen(vm["outfile"].as<string>().c_str(), "w");
        if (!outfile) throw std::ios_base::failure(strerror(errno), std::error_code(errno, std::system_category()));
        outputs.push_back(outfile);
    }
    FILE* statFile = nullptr;
    if (vm.count("stat-file")) {
        auto statFileName = vm["stat-file"].as<string>();
        fs::create_directories(fs::absolute(fs::path(statFileName)).parent_path());
        statFile = fopen(statFileName.c_str(), "w");
        if (!statFile) throw std::ios_base::failure(strerror(errno), std::error_code(errno, std::system_category()));
        fmt::print(statFile, "#{}\n", fmt::join(std::span(argv, argc + argv), " "));
        fmt::print(statFile, "# degrees format: <#zi+inactive>:<#zi+blank>:<#zi+active>:<#nonzi+inactive>:<#nonzi+blank>:<#nonzi.active>;...\n");
        fmt::print(statFile, "{}", "name,n,m,n_zero_injection,n_pmu,n_inactive,n_blank,n_observed,propagation_distance,tree_width,degrees\n");
    }
    FILE* fortStatFile = nullptr;
    if (vm.count("fort-stats")) {
        auto fortStatFileName = vm["fort-stats"].as<string>();
        fs::create_directories(fs::absolute(fs::path(fortStatFileName)).parent_path());
        fortStatFile = fopen(fortStatFileName.c_str(), "w");
        if (!fortStatFile) throw std::ios_base::failure(strerror(errno), std::error_code(errno, std::system_category()));
        fmt::print(fortStatFile, "#{}\n", fmt::join(std::span(argv, argc + argv), " "));
        fmt::print(fortStatFile, "{}", "name,run,forts,avg_size,num_hitting_sets\n");
    }

    for (auto out: outputs) {
        fmt::print(out, "#{}\n", fmt::join(std::span(argv, argc + argv), " "));
        fmt::print(out, "{}\n","name,run,bound,pmus,solved,result,t_total,t_reductions,t_solver,n,m,zi,n_reduced,m_reduced,zi_reduced,pmu_reduced,blank_reduced");
    }

    for (const std::string& filename: inputs) {
        string currentName = fs::path(filename).filename().string();
        currentName = currentName.substr(0, currentName.rfind('.'));
        PdsState inputState(readAutoGraph(filename, allZeroInjection));
        for (size_t run = 0; run < repetitions; ++run) {
            struct BoundInfo {
                size_t lower, upper;
                size_t extraInfo;
                long time;
            };
            std::vector<BoundInfo> bounds;
            auto state = inputState;
            size_t n = state.graph().numVertices();
            size_t m = state.graph().numEdges();
            size_t zi = state.numZeroInjection();
            if (drawdir) {
                writePds(state.graph(), fmt::format("{}/0_input.pds", drawdir->string()));
            }

            auto t0 = now();
            reduce(state);
            size_t nReduced = state.graph().numVertices();
            size_t mReduced = state.graph().numEdges();
            size_t ziReduced = state.numZeroInjection();
            size_t pmusReduced = state.numActive();
            size_t blankReduced = nReduced - state.numActive() - state.numInactive();
            SolveResult result = {state.numActive(), state.numActive() + state.numBlank(), SolveState::Optimal };
            auto reduced = state;
            FILE* boundsFile = nullptr;
            if (vm.count("write-bounds")) {
                std::string boundsDir = vm["write-bounds"].as<string>();
                fs::create_directories(fs::absolute(fs::path(boundsDir)));
                auto boundsFileName = fmt::format("{}/{}-{}-bounds.csv", boundsDir, currentName, run);
                boundsFile = fopen(boundsFileName.c_str(), "w+");
                if (!boundsFile) throw std::ios_base::failure(strerror(errno), std::error_code(errno, std::system_category()));
                fmt::print(boundsFile, "#{}\n", fmt::join(std::span(argv, argc + argv), " "));
                fmt::print(boundsFile, "name,run,subinstance,t,lower,upper,gap,extra\n");
            }
            auto t1 = now();
            if (state.allObserved()) {
                if (boundsFile) {
                    const auto &b = BoundInfo{state.numActive(), state.numActive(), 0, µs(now() - t0)};
                    double gap = (double(b.upper) - double(b.lower)) / double(b.upper);
                    fmt::print(boundsFile, "{},{},0,{},{},{},{},0\n", currentName, run, b.time, b.lower, b.upper, gap);
                    fflush(boundsFile);
                }
                result = SolveResult { state.numActive(), state.numActive(), SolveState::Optimal};
            } else if (subproblems) {
                auto checkpoint = t1;
                auto subproblems = state.subproblems();
                ranges::sort(subproblems, [](const auto& left, const auto& right) { return left.graph().numVertices() < right.graph().numVertices(); });
                size_t subproblemNumber = 0;
                for (auto& substate: subproblems) {
                    if (!substate.allObserved()) {
                        size_t initialActive = substate.numActive();
                        size_t initialBlank = substate.numBlank();
                        result.lower -= initialActive;
                        result.upper -= initialBlank + initialActive;
                        auto boundCB = [&bounds, t0, result, boundsFile, &currentName, &subproblemNumber, run, initialBlank, initialActive](size_t lower, size_t upper, size_t extra) {
                            auto time = now();
                            bounds.emplace_back(BoundInfo{
                                result.lower + std::max(lower, initialActive),
                                result.upper + std::min(upper, initialBlank + initialActive),
                                extra, µs(time - t0)
                            });
                            if (boundsFile) {
                                const auto& b = bounds.back();
                                double gap = (double(b.upper) - double(b.lower)) / double(b.upper);
                                fmt::print(boundsFile, "{},{},{},{},{},{},{},{}\n", currentName, run, subproblemNumber, b.time, b.lower, b.upper, gap, b.extraInfo);
                                fflush(boundsFile);
                            }
                        };
                        Solver solve = getSolver(vm, fortStats, cycleStats, pathStats, currentName, run, subproblemNumber, boundCB);
                        auto tSubproblem = now();
                        double remainingTimeout = std::max(1.0, timeout - std::chrono::duration_cast<std::chrono::seconds>(tSubproblem - checkpoint).count());
                        auto subresult = solve(substate, remainingTimeout);
                        result.state = combineSolveState(result.state, subresult.state);
                        state.applySubsolution(substate);
                        result.lower += std::max(subresult.lower, initialActive);
                        result.upper += std::max(subresult.upper, initialActive);
                        boundCB(subresult.lower, subresult.upper, 0);
                    }
                    ++subproblemNumber;
                }
            } else {
                auto boundCB = [&bounds, t0, boundsFile, &currentName, run](size_t lower, size_t upper, size_t extra) {
                    auto time = now();
                    bounds.emplace_back(BoundInfo{lower, upper, extra, µs(time - t0)});
                    if (boundsFile) {
                        const auto& b = bounds.back();
                        double gap = (double(b.upper) - double(b.lower)) / double(b.upper);
                        fmt::print(boundsFile, "{},{},{},{},{},{},{},{}\n", currentName, run, 0, b.time, b.lower, b.upper, gap, b.extraInfo);
                        fflush(boundsFile);
                    }
                };
                Solver solve = getSolver(vm, fortStats, cycleStats, pathStats, currentName, run, 0, boundCB);
                result = solve(state, timeout);
                boundCB(result.lower, result.upper, -1);
            }
            auto t2 = now();
            if (drawdir) {
                writePds(reduced.graph(), fmt::format("{}/1_reductions.pds", drawdir->string()));
                writePds(state.graph(), fmt::format("{}/2_solved_preprocessed.pds", drawdir->string()));
            }
            if (solDir) {
                auto name = fs::path(filename).filename().string();
                auto end = name.rfind('.');
                name = name.substr(0, end);
                auto solPath = *solDir / fmt::format("{}_{}.pds", name, run);
                auto solution = inputState.graph();
                for (auto v: solution.vertices()) {
                    if (!state.graph().hasVertex(v)) {
                        solution.getVertex(v).pmu = PmuState::Inactive;
                    } else {
                        if (state.isActive(v)) {
                            solution.getVertex(v).pmu = PmuState::Active;
                        } else if (state.isInactive(v)) {
                            solution.getVertex(v).pmu = PmuState::Inactive;
                        }
                    }
                }
                writePds(solution, solPath);
            }
            size_t pmus = state.numActive();
            fmt::memory_buffer buf;
            using namespace fmt::literals;
            for (auto file: outputs) {
                fmt::print(file,
                           "{name},{run},{lower_bound},{pmus},{solved},{result},{t_total},{t_reductions},{t_solver},{n},{m},{zi},{nReduced},{mReduced},{ziReduced},{pmuReduced},{blankReduced}\n",
                           "name"_a = filename,
                           "run"_a = run,
                           "lower_bound"_a = result.lower,
                           "pmus"_a = pmus,
                           "solved"_a = state.allObserved(),
                           "result"_a = result.state,
                           "t_total"_a = µs(t2 - t0),
                           "t_reductions"_a = µs(t1 - t0),
                           "t_solver"_a = µs(t2 - t1),
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
            if (statFile != nullptr) {
                writeSolutionStatistics(filename, state, statFile);
            }
            if (fortStatFile != nullptr) {
                writeFortStats(filename, run, fortStats, fortStatFile);
                writeFortStats(filename, run, cycleStats, fortStatFile);
                writeFortStats(filename, run, pathStats, fortStatFile);
            }
        }
    }
    if (statFile != nullptr) {
        fclose(statFile);
    }
    if (outfile != nullptr) {
        fclose(outfile);
    }
}