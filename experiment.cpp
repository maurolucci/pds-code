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

int main(int argc, const char** argv) {
    namespace po = boost::program_options;
    namespace fs = std::filesystem;
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
            ("repeat,n", po::value<size_t>()->default_value(1)->implicit_value(5), "number of experiment repetitions")
            ("solve,s", po::value<string>()->default_value("none")->implicit_value("gurobi"), "solve method. can be any of [none,gurobi,greedy]")
            ("subproblems,u", "split into subproblems before calling solve")
            ("timeout,t", po::value<double>()->default_value(600.), "gurobi time limit (seconds)")
            ("draw,d", po::value<string>()->implicit_value("out"s), "draw states")
            ("write,w", po::value<string>()->implicit_value("solutions"s), "write solutions to the specified directory")
            ("stat-file", po::value<string>(), "write statistics about solutions")
            ("greedy-bounds,b", po::value<int>()->default_value(0)->implicit_value(1), "if possible, use a greedy algorithm to compute an upper bound (0: never, 1: without reductions (faster), 2: with reductions (more precise))")
            ("fort-stats", po::value<string>(), "file for fort statistics")
            ("early-stop", "stop hitting set solver when violating hitting set is found")
            ("write-forts", po::value<string>()->implicit_value("hs"), "directory to which to write hitting set instance")
            ("verbose,v", "print additional solver status info")
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

    string solverName = vm["solve"].as<string>();
    std::function<SolveResult(PdsState&, double)> solve;
    bool verbose = vm.count("verbose");
    bool earlyStop = vm.count("early-stop");
    preloadMIPSolver();
    int greedyBoundMode = vm["greedy-bounds"].as<int>();
    size_t subproblemNumber = 0;
    string currentName;
    if (vm.count("write-forts")) {
        std::string fortsDirName = vm["write-forts"].as<string>();
        fs::create_directories(fs::absolute(fs::path(fortsDirName)));
    }
    FortStats fortStats;
    callback::FortCallback fortCallback = [&,solved=size_t{0}] (callback::When when, const PdsState& state, const std::vector<VertexList>& forts, size_t lower, size_t upper) mutable {
        if (vm.count("fort-stats") && when == pds::callback::When::FINAL) {
            size_t totalFortSize = 0;
            for (auto& fort: forts) {
                totalFortSize += fort.size();
            }
            double averageSize = double(totalFortSize) / double(forts.size());
            fortStats = {averageSize, forts.size(), solved};
        }
        if (vm.count("write-forts")) {
            {
                auto hs_file = fmt::format("hs/{}-{}-hs-{:04}.hs", currentName, subproblemNumber, solved);
                FILE* file = fopen(hs_file.c_str(), "w");
                fmt::print(file, "{} {}\n", state.graph().numVertices(), forts.size() + state.numActive());
                // ensure that active vertices are selected to get the correct bound
                for (auto v: state.graph().vertices()) {
                    if (state.isActive(v)) {
                        fmt::print(file, "1 {}\n", v);
                    }
                }
                for (auto& f: forts) {
                    VertexList fortBlank;
                    bool skip = false;
                    for (auto v: f) {
                        if (!state.isInactive(v)) {
                            fortBlank.push_back(v);
                        }
                    }
                    if (!skip) {
                        ranges::sort(fortBlank);
                        fmt::print(file, "{} {}\n", fortBlank.size(), fmt::join(fortBlank, " "));
                    }
                }
                fclose(file);
            }
        }
        ++solved;
    };
    if (solverName == "branching") {
        solve = [](auto& state, double) { return solveBranching(state, true, greedy_strategies::largestDegree); };
    } else if (solverName == "greedy") {
        solve = [](auto& state, double) { return solveGreedy(state, true, greedy_strategies::largestDegree); };
    } else if (solverName == "fast") {
        solve = [](auto& state, double) { return fastGreedy(state, true); };
    } else if (solverName == "topdown") {
        solve = [](auto& state, double) { return topDownGreedy(state); };
    } else if (solverName == "none") {
        solve = [](auto & state, double) { return SolveResult{ size_t{0}, state.numActive() + state.numBlank(), SolveState::Other }; };
    } else if (solverName == "smith") {
        solve = [=](auto& state, double timeLimit) {
            return solveBozeman(state, verbose, timeLimit, 1, greedyBoundMode, earlyStop, fortCallback);
        };
    } else if (solverName == "bozeman") {
        solve = [=](auto& state, double timeLimit) {
            return solveBozeman(state, verbose, timeLimit, 0, greedyBoundMode, earlyStop, fortCallback);
        };
    } else if (solverName == "bozeman2") {
        solve = [=](auto& state, double timeLimit) {
            return solveBozeman(state, verbose, timeLimit, 2, greedyBoundMode, earlyStop, fortCallback);
        };
    } else if (solverName == "bozeman3") {
        solve = [=](auto& state, double timeLimit) {
            return solveBozeman(state, verbose, timeLimit, 3, greedyBoundMode, earlyStop, fortCallback);
        };
    } else if (solverName == "forts") {
        solve = [=](auto& state, double timeLimit) {
            return solveBozeman(state, verbose, timeLimit, 4, greedyBoundMode, earlyStop, fortCallback);
        };
    } else {
        try {
            solve = [model=getModel(solverName),verbose](auto &state, double timeout) {
                return solvePowerDominatingSet(state, verbose, timeout, model);
            };
        } catch(std::invalid_argument& ex) {
            fmt::print(stderr, "{}", ex.what());
            return 2;
        }
    }

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
        outputs.push_back(outfile);
    }
    FILE* statFile = nullptr;
    if (vm.count("stat-file")) {
        auto statFileName = vm["stat-file"].as<string>();
        fs::create_directories(fs::absolute(fs::path(statFileName)).parent_path());
        statFile = fopen(statFileName.c_str(), "w");
        fmt::print(statFile, "#{}\n", fmt::join(std::span(argv, argc + argv), " "));
        fmt::print(statFile, "# degrees format: <#zi+inactive>:<#zi+blank>:<#zi+active>:<#nonzi+inactive>:<#nonzi+blank>:<#nonzi.active>;...\n");
        fmt::print(statFile, "{}", "name,n,m,n_zero_injection,n_pmu,n_inactive,n_blank,n_observed,propagation_distance,tree_width,degrees\n");
    }
    FILE* fortStatFile = nullptr;
    if (vm.count("fort-stats")) {
        auto fortStatFileName = vm["fort-stats"].as<string>();
        fs::create_directories(fs::absolute(fs::path(fortStatFileName)).parent_path());
        fortStatFile = fopen(fortStatFileName.c_str(), "w");
        fmt::print(fortStatFile, "#{}\n", fmt::join(std::span(argv, argc + argv), " "));
        fmt::print(fortStatFile, "{}", "name,run,forts,avg_size,num_hitting_sets\n");
    }

    for (auto out: outputs) {
        fmt::print(out, "#{}\n", fmt::join(std::span(argv, argc + argv), " "));
        fmt::print(out, "{}\n","name,run,pmus,solved,result,t_total,t_reductions,t_solver,n,m,zi,n_reduced,m_reduced,zi_reduced,pmu_reduced,blank_reduced");
    }

    for (const std::string& filename: inputs) {
        currentName = fs::path(filename).filename().string();
        currentName = currentName.substr(0, currentName.rfind('.'));
        PdsState inputState(readAutoGraph(filename, allZeroInjection));
        for (size_t run = 0; run < repetitions; ++run) {
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
            auto t1 = now();
            SolveResult result = {state.numActive(), state.numActive(), SolveState::Optimal };
            auto reduced = state;
            if (subproblems) {
                auto checkpoint = t1;
                auto subproblems = state.subproblems();
                ranges::sort(subproblems, [](const auto& left, const auto& right) { return left.graph().numVertices() < right.graph().numVertices(); });
                for (auto& substate: subproblems) {
                    if (!substate.allObserved()) {
                        auto tSubproblem = now();
                        double remainingTimeout = std::max(1.0, timeout - std::chrono::duration_cast<std::chrono::seconds>(tSubproblem - checkpoint).count());
                        size_t initialActive = substate.numActive();
                        auto subresult = solve(substate, remainingTimeout);
                        result.state = combineSolveState(result.state, subresult.state);
                        state.applySubsolution(substate);
                        result.lower += std::max(subresult.lower, initialActive) - initialActive;
                        result.upper += std::max(subresult.upper, initialActive) - initialActive;
                    }
                }
            } else {
                result = solve(state, timeout);
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