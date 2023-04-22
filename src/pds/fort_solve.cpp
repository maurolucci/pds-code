#include "fort_solve.hpp"

#include <gurobi_c++.h>
#include <range/v3/all.hpp>

#include "pdssolve.hpp"
#include "gurobi_solve.hpp"

namespace pds {


VertexList bozemanFortNeighborhood(const PdsState &state, bool output, double timeLimit, VertexSet& seen) {
    GRBModel model(getEnv());
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    unused(state, output, timeLimit);
    VertexMap<GRBVar> xi;
    for (auto v: state.graph().vertices()) {
        xi.emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("f_{}", v)));
    }
    //model.setObjective(GRBLinExpr{0.0});
    GRBLinExpr allXi;
    for (auto v: state.graph().vertices()) {
        if (state.isObserved(v)) {
            model.addConstr(xi[v] == 0);
        }
        allXi += xi[v];
        for (auto w: state.graph().neighbors(v)) {
            if (!state.isZeroInjection(w)) continue;
            GRBLinExpr sum;
            for (auto u: state.graph().neighbors(w)) {
                if (u != v) {
                    sum += xi[u];
                }
            }
            model.addConstr(xi[w] + sum >= xi[v]);
        }
    }
    model.addConstr(allXi >= 1);
    model.optimize();
    VertexList fort;
    switch (model.get(GRB_IntAttr_Status)) {
        case GRB_OPTIMAL:
        case GRB_TIME_LIMIT:
            assert(seen.empty());
            for (auto v: state.graph().vertices()) {
                if (xi[v].get(GRB_DoubleAttr_X) > 0.5) {
                    if (!seen.contains(v)) fort.push_back(v);
                    for (auto w: state.graph().neighbors(v)) {
                        if (!seen.contains(w)) fort.push_back(w);
                    }
                }
            }
            for (auto v: fort) seen.erase(v);
            assert(seen.empty());
            return fort;
        default:
            return {};
    }
}
VertexList bozemanFortNeighborhood2(const PdsState &state, bool output, double timeLimit) {
    GRBModel model(getEnv());
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    //model.set(GRB_StringParam_LogFile, "gurobi.log");
    //model.set(GRB_DoubleParam_MIPGap, 0.05);
    //model.set(GRB_DoubleParam_MIPGapAbs, 10);
    unused(state, output, timeLimit);
    VertexMap<GRBVar> xi;
    VertexMap<GRBVar> ni;
    for (auto v: state.graph().vertices()) {
        xi.emplace(v, model.addVar(0.0, 1.0, 0.0, GRB_BINARY, fmt::format("f_{}", v)));
        ni.emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("n_{}", v)));
    }
    //model.setObjective(GRBLinExpr{0.0});
    GRBLinExpr allXi;
    for (auto v: state.graph().vertices()) {
        if (state.isObserved(v)) {
            model.addConstr(xi[v] == 0);
        }
        allXi += xi[v];
        GRBLinExpr neighborSum{xi.at(v)};
        for (auto w: state.graph().neighbors(v)) {
            neighborSum += xi.at(w);
            if (!state.isZeroInjection(w)) continue;
            GRBLinExpr sum;
            for (auto u: state.graph().neighbors(w)) {
                if (u != v) {
                    sum += xi[u];
                }
            }
            model.addConstr(xi[w] + sum >= xi[v]);
            //model.addConstr(ni.at(v) >= xi.at(w));
        }
        //model.addConstr(ni.at(v) >= xi.at(v));
        model.addConstr(neighborSum <= state.graph().numVertices() * ni.at(v));
    }
    model.addConstr(allXi >= 1);
    model.optimize();
    VertexList fort;
    switch (model.get(GRB_IntAttr_Status)) {
        case GRB_OPTIMAL:
        case GRB_TIME_LIMIT:
            for (auto v: state.graph().vertices()) {
                if (ni.at(v).get(GRB_DoubleAttr_X) > 0.5) {
                    fort.push_back(v);
                }
            }
            return fort;
        default:
            return {};
    }
}

VertexList bozemanFortNeighborhood3(const PdsState &state, bool output, double timeLimit) {
    GRBModel model(getEnv());
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    //model.set(GRB_DoubleParam_MIPGap, 0.05);
    //model.set(GRB_DoubleParam_MIPGapAbs, 10);
    unused(state, output, timeLimit);
    VertexMap<GRBVar> xi;
    VertexMap<GRBVar> ni;
    for (auto v: state.graph().vertices()) {
        double weight = state.isInactive(v) ? 0.0 : 1.0;
        xi.emplace(v, model.addVar(0.0, 1.0, 0.0, GRB_BINARY, fmt::format("f_{}", v)));
        ni.emplace(v, model.addVar(0.0, 1.0, weight, GRB_BINARY, fmt::format("n_{}", v)));
    }
    //model.setObjective(GRBLinExpr{0.0});
    GRBLinExpr allXi;
    for (auto v: state.graph().vertices()) {
        if (state.isObserved(v)) {
            model.addConstr(xi[v] == 0);
        }
        if (state.isBlank(v)) {
            allXi += xi[v];
        }
        GRBLinExpr neighborSum{xi.at(v)};
        for (auto w: state.graph().neighbors(v)) {
            neighborSum += xi.at(w);
            if (!state.isZeroInjection(w)) continue;
            GRBLinExpr sum;
            for (auto u: state.graph().neighbors(w)) {
                if (u != v) {
                    sum += xi[u];
                }
            }
            model.addConstr(xi[w] + sum >= xi[v]);
        }
        model.addConstr(neighborSum <= state.graph().numVertices() * ni.at(v));
    }
    model.addConstr(allXi >= 1);
    model.optimize();
    VertexList fort;
    switch (model.get(GRB_IntAttr_Status)) {
        case GRB_OPTIMAL:
        case GRB_TIME_LIMIT:
            for (auto v: state.graph().vertices()) {
                if (ni.at(v).get(GRB_DoubleAttr_X) > 0.5) {
                    fort.push_back(v);
                }
            }
            return fort;
        default:
            return {};
    }
}

namespace { struct Unimplemented {}; }
std::vector<VertexList> components(const PdsState& state, VertexSet& seen) {
    const auto& graph = state.graph();
    std::vector<VertexList> components;
    VertexList stack;
    for (auto start: graph.vertices()) {
        if (!seen.contains(start)) {
            stack.push_back(start);
            seen.insert(start);
            VertexList comp;
            while (!stack.empty()) {
                auto current = stack.back(); stack.pop_back();
                comp.push_back(current);
                for (auto w: graph.neighbors(current)) {
                    if (!seen.contains(w)) {
                        seen.insert(w);
                        stack.push_back(w);
                    }
                }
            }
            components.push_back(comp);
        }
    }
    for (auto& comp: components) {
        for (auto v: comp) {
            seen.erase(v);
        }
    }
    return components;
}

namespace {
auto now() {
    return std::chrono::high_resolution_clock::now();
}
}

VertexList loganFortNeighborhood(const PdsState& state, bool output, double timeLimit, VertexSet& seen) {
    VertexSet junctions;
    std::vector<PowerGrid::VertexDescriptor> junctionVertices;
    for (auto v: state.graph().vertices()) {
        if (state.graph().degree(v) + state.isZeroInjection(v) >= 3) {
            junctions.insert(v);
            seen.insert(v);
            junctionVertices.push_back(v);
        }
    }
    auto comp = components(state, seen);
    for (auto v: junctionVertices) {
        seen.erase(v);
    }
    GRBModel model(getEnv());
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    VertexMap<GRBVar> fv;
    VertexMap<GRBVar> mv;
    std::vector<GRBVar> fp;
    for (auto v: junctionVertices) {
        fv.emplace(v, model.addVar(0.0, 1.0, 0.0, GRB_BINARY, fmt::format("F_{}", v)));
        mv.emplace(v, model.addVar(0.0, 1.0, 0.0, GRB_BINARY, fmt::format("M_{}", v)));
    }
    for (size_t i = 0; i < comp.size(); ++i) {
        fp.push_back(model.addVar(0.0, 1.0, 0.0, GRB_BINARY, fmt::format("F_p{}", i)));
    }
    GRBLinExpr sumMvFp;
    GRBLinExpr objective;
    GRBLinExpr weightedSumMvFp;
    GRBLinExpr activeSum;
    for (auto v: junctionVertices) {
        sumMvFp += mv.at(v);
        if (state.isObserved(v)) {
            weightedSumMvFp += mv.at(v);
        }
        if (state.isActive(v)) {
            activeSum += mv.at(v);
        }
    }
    for (size_t i = 0; i < comp.size(); ++i) {
        sumMvFp += comp.size() * fp.at(i);
        for (auto v: comp[i]) {
            if (state.isObserved(v)) { weightedSumMvFp += fp.at(i); }
            if (state.isActive(v)) { activeSum += fp.at(i); }
        }
    }
    model.setObjective(sumMvFp);
    //model.addConstr(-weightedSumMvFp == 0); // 7
    model.addConstr(activeSum == 0);
    model.addConstr(sumMvFp >= 1);
    for (auto u: junctionVertices) {
        model.addConstr(fv.at(u) <= mv.at(u));
        for (auto v: state.graph().neighbors(u)) {
            if (junctions.contains(v)) {
                model.addConstr(fv.at(v) <= mv.at(u));
            }
        }
    }
    for (size_t i = 0; i < comp.size(); ++i) {
        const auto& path = comp[i];
        for (auto u: path) {
            for (auto v: state.graph().neighbors(u)) {
                if (seen.contains(v)) continue;
                seen.insert(v);
                if (junctions.contains(v)) {
                    model.addConstr(fv.at(v) <= fp.at(i));
                    model.addConstr(fp.at(i) <= mv.at(v));
                }
            }
        }
        for (auto v: path) {
            seen.erase(v);
            for (auto w: state.graph().neighbors(v)) {
                seen.erase(w);
            }
        }
    }
    assert(seen.empty());
    for (auto v: junctionVertices) {
        GRBLinExpr neighborSum;
        for (auto u: state.graph().neighbors(v)) {
            if (junctions.contains(u)) { neighborSum += fv.at(u); }
        }
        for (size_t i = 0; i < comp.size(); ++i) {
            const auto& path = comp[i];
            for (auto p: path) { seen.insert(p); }
            size_t common = 0;
            for (auto w: state.graph().neighbors(v)) {
                if (seen.contains(w)) { ++common; }
            }
            if (common) { neighborSum += common * fp.at(i); }
            for (auto p: path) { seen.erase(p); }
        }
        model.addConstr(2 * (mv.at(v) - fv.at(v)) <= neighborSum);
    }
    assert(seen.empty());
    model.optimize();
    VertexList fort;
    switch (model.get(GRB_IntAttr_Status)) {
        case GRB_OPTIMAL:
        case GRB_TIME_LIMIT:
            for (auto v: junctionVertices) {
                if (mv.at(v).get(GRB_DoubleAttr_X) > 0.5) {
                    fort.push_back(v);
                }
            }
            for (size_t i = 0; i < comp.size(); ++i) {
                if (fp.at(i).get(GRB_DoubleAttr_X) > 0.5) {
                    for (auto v: comp[i]) { fort.push_back(v); }
                }
            }
            return fort;
        default:
            return {};
    }
}

std::vector<VertexList> initialForts(const PdsState& state, VertexSet& seen) {
    std::vector<VertexList> forts;
    assert(seen.empty());
    for (auto v: state.graph().vertices()) { seen.insert(v); }
    for (auto v: state.graph().vertices()) {
        if (!state.isObserved(v)) {
            seen.erase(v);
            for (auto w: state.graph().neighbors(v)) {
                seen.erase(w);
            }
        }
    }
    auto unobserved = components(state, seen);
    for (auto v: state.graph().vertices()) { seen.erase(v); }
    assert(seen.empty());
    for (auto& comp: unobserved) {
        VertexList fort;
        assert(seen.empty());
        for (auto v: comp) {
            if (!seen.contains(v)) {
                fort.push_back(v);
                seen.insert(v);
            }
        }
        for (auto v: fort) {
            seen.erase(v);
        }
        assert(seen.empty());
        forts.emplace_back(std::move(fort));
    }
    return forts;
}

std::vector<VertexList> initialForts2(PdsState& state, VertexSet& seen) {
    auto blank = state.graph().vertices()
                 | ranges::views::filter([&state](auto v) { return state.isBlank(v); })
                 | ranges::to<std::vector>;
    //ranges::shuffle(blank);
    for (auto v: blank) { state.setActive(v); }
    size_t first_blank = 0;
    std::vector<VertexList> forts;
    for (size_t i = 0; i < blank.size(); ++i) {
        assert(state.isActive(blank[i]));
        state.setBlank(blank[i]);
        if (!state.allObserved()) {
            VertexList fort;
            for (auto v: state.graph().neighbors(blank[i])) {
                if (!state.isObserved(v)) {
                    seen.insert(v);
                    fort.push_back(v);
                }
            }
            if (!state.isObserved(blank[i])) {
                seen.insert(blank[i]);
                fort.push_back(blank[i]);
            }
            assert(!fort.empty());
            size_t start = 0;
            while (start < fort.size()) {
                auto current = fort[start]; ++start;
                for (auto other: state.graph().neighbors(current)) {
                    if (!seen.contains(other) && (!state.isObserved(current) || !state.isObserved(other))) {
                        seen.insert(other);
                        fort.push_back(other);
                    }
                }
            }
            for (auto v: fort) { seen.erase(v); }
            for (; first_blank <= i; ++first_blank) {
                state.setActive(blank[first_blank]);
            }
            forts.emplace_back(std::move(fort));
        }
    }
    for (auto v: blank) { state.setBlank(v); }
    return forts;
}

template<class Span>
size_t localSearchFortShrinker(PdsState& state, Span blank) {
    constexpr size_t INVALID = -1;
    size_t currentMin = state.numUnobserved();
    size_t blank_size = blank.size();
    size_t bestVertex = INVALID;
    do {
        bestVertex = INVALID;
        currentMin = state.numUnobserved();
        for (size_t i = blank_size; i--;) {
            state.setActive(blank[i]);
            if (!state.allObserved()) {
                if (state.numUnobserved() <= currentMin) {
                    bestVertex = i;
                }
            }
            state.setBlank(blank[i]);
        }
        if (bestVertex != INVALID) {
            state.setActive(blank[bestVertex]);
            --blank_size;
            std::swap(blank[bestVertex], blank[blank_size]);
        }
    } while (bestVertex != INVALID);
    return blank.size() - blank_size;
}

template<class Span>
size_t greedyFortShrinker(PdsState& state, Span deselected) {
    if (deselected.empty()) return 0;
    ranges::shuffle(deselected);
    size_t skip = 0;
    size_t vertex_count = deselected.size();
    for (size_t i = vertex_count; i--; ) {
        state.setActive(deselected[i]);
        if (!state.allObserved()) {
            ++skip;
            std::swap(deselected[deselected.size() - skip], deselected[i]);
        } else {
            state.setBlank(deselected[i]);
        }
    }
    return skip;
}

inline VertexList findFort(PdsState& state, PdsState::Vertex start, VertexSet& seen) {
    VertexList fort;
    seen.insert(start);
    fort.push_back(start);
    assert(!fort.empty());
    size_t front = 0;
    while (front < fort.size()) {
        auto current = fort[front];
        ++front;
        for (auto other: state.graph().neighbors(current)) {
            if (!seen.contains(other) && (!state.isObserved(current) || !state.isObserved(other))) {
                seen.insert(other);
                fort.push_back(other);
            }
        }
    }
    for (auto v: fort) { seen.erase(v); }
    return fort;
}

std::vector<VertexList> initialForts3(PdsState& state, VertexSet& seen) {
    auto blank = state.graph().vertices()
                 | ranges::views::filter([&state](auto v) { return state.isBlank(v); })
                 | ranges::to<std::vector>;
    for (auto v: blank) { state.setActive(v); }
    std::vector<VertexList> forts;
    ranges::shuffle(blank);
    size_t blank_size = blank.size();
    size_t i = 0;
    while (i < blank_size) {
        state.setBlank(blank[i]);
        if (!state.allObserved()) {
            auto fort = findFort(state, blank[i], seen);
            assert(!fort.empty());
            assert(fort.front() == blank[i]);
            size_t reselected = localSearchFortShrinker(state, std::span(fort.begin(), fort.size()));
            if (reselected) {
                forts.emplace_back(findFort(state, blank[i], seen));
                for (size_t j = 1; j <= reselected; ++j) {
                    state.setBlank(fort[fort.size() - j]);
                }
            } else {
                forts.emplace_back(std::move(fort));
            }
            state.setActive(blank[i]);
            --blank_size;
            std::swap(blank[i], blank[blank_size]);
        } else {
            ++i;
        }
    }
    for (auto v: blank) { state.setBlank(v); }
    return forts;
}

namespace {
struct Callback : public GRBCallback {
    size_t* lower;
    size_t* upper;
    VertexMap<GRBVar>* pi;
    const PdsState* base;
    PdsState* solution;
    PdsState* upperBound;
    std::span<PdsState::Vertex> blank;
    bool earlyStop;
    Callback(size_t& lower, size_t& upper, VertexMap<GRBVar>& pi,
             const PdsState& base, PdsState& solution, PdsState & upperBound,
             std::span<PdsState::Vertex> blank, bool earlyStop)
            : lower(&lower), upper(&upper), pi(&pi), base(&base), solution(&solution),
              upperBound(&upperBound), blank(blank), earlyStop(earlyStop)
    { }
    void callback() override {
        if (where == GRB_CB_MIP) {
            if (getIntInfo(GRB_CB_MIP_SOLCNT) > 0) {
                auto objVal = static_cast<size_t>(getDoubleInfo(GRB_CB_MIP_OBJBST) + 0.5);
                auto objBound = static_cast<size_t>(getDoubleInfo(GRB_CB_MIP_OBJBND));
                if (objVal <= *upper && !solution->allObserved() && earlyStop && objBound > *lower) {
                    abort();
                }
            }
        }
        if (where == GRB_CB_MIPSOL) {
            auto objVal = static_cast<size_t>(getDoubleInfo(GRB_CB_MIPSOL_OBJBST) + 0.5);
            auto objBound = static_cast<size_t>(getDoubleInfo(GRB_CB_MIPSOL_OBJBND));
            if (objVal <= *upper) {
                for (auto v: blank) {
                    if (getSolution(pi->at(v)) > 0.5) {
                        solution->setActive(v);
                    } else if (base->isBlank(v)) {
                        solution->setBlank(v);
                    } else {
                        solution->setInactive(v);
                    }
                }
                //fmt::print("feasible solution {} <= {}; {} <= {}; {}\n", *lower, objBound, objVal, *upper, solution->allObserved());
                if (!solution->allObserved()) {
                    if (earlyStop && objBound > *lower) { abort(); }
                } else if (objVal < *upper) {
                    fmt::print("gurobi H {}\n", objVal);
                    *upper = objVal;
                    *upperBound = *solution;
                    if (*upper <= *lower) { abort(); }
                }
            }
        }
    }
};
}

SolveResult solveBozeman(
        PdsState &state,
        int output,
        double timeLimit,
        int variant,
        int greedyUpper,
        bool earlyStop,
        callback::FortCallback callback
) {
    unused(state, output, timeLimit);
    auto lastSolution = state;
    //writePds(lastSolution.graph(), fmt::format("out/0_input.pds"));
    auto blankVertices = state.graph().vertices()
                         | ranges::views::filter([&state] (auto v) { return state.isBlank(v); })
                         | ranges::to<std::vector<PdsState::Vertex>>;
    PdsState feasibleSolution = state;

    size_t lowerBound = 0;
    size_t upperBound = state.numActive() + state.numBlank();
    try {
        VertexSet seen;
        auto forts = initialForts2(state, seen);
        GRBModel model(getEnv());
        model.set(GRB_IntParam_LogToConsole, false);
        model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        model.set(GRB_IntParam_NumericFocus, 1);
        model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
        if (output) {
            model.set(GRB_StringParam_LogFile, "gurobi.log");
        }
        VertexMap<GRBVar> pi;
        for (auto v: state.graph().vertices()) {
            pi.emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("p_{}", v)));
            if (state.isActive(v)) {
                model.addConstr(pi.at(v) == 1);
            }
            if (state.isInactive(v)) {
                model.addConstr(pi.at(v) == 0);
            }
        }
        Callback cb(lowerBound, upperBound, pi, state, lastSolution, feasibleSolution, blankVertices, earlyStop);
        model.setCallback(&cb);
        auto startingTime = now();
        int status;
        size_t processedForts = 0;
        size_t totalFortSize = 0;
        while (true) {
            auto currentTime = now();
            double remainingTimeout = std::max(1.0, timeLimit - std::chrono::duration_cast<std::chrono::seconds>(currentTime - startingTime).count());
            switch(variant) {
                case 1:
                    forts.emplace_back(loganFortNeighborhood(lastSolution, false, remainingTimeout, seen));
                    break;
                case 2:
                    forts.emplace_back(bozemanFortNeighborhood2(lastSolution, false, remainingTimeout));
                    break;
                case 3:
                    forts.emplace_back(bozemanFortNeighborhood3(lastSolution, false, remainingTimeout));
                    break;
                case 4: {
                    if (output) { fmt::print("finding forts in {} blank vertices\n", lastSolution.numBlank()); }
                    auto more_forts = initialForts3(lastSolution, seen);
                    for (auto f: more_forts) {
                        forts.emplace_back(std::move(f));
                    }
                    break; }
                case 0:
                default:
                    forts.emplace_back(bozemanFortNeighborhood(lastSolution, output, remainingTimeout, seen));
            }
            if (forts.back().empty()) break;
            for (; processedForts < forts.size(); ++processedForts) {
                GRBLinExpr fortSum;
                size_t blank = 0;
                for (auto& v: forts[processedForts]) {
                    if (lastSolution.isActive(v)) {
                        if (output) {
                            fmt::print("!!!active vertex in neighborhood!!! {} {} {}\n",
                                       v,
                                       lastSolution.numObserved(),
                                       lastSolution.graph().numVertices());
                        }
                    }
                    if (state.isBlank(v)) {
                        fortSum += pi.at(v);
                        ++blank;
                    } else if (state.isInactive(v)) {
                        // might have been changed by initialFortNeighborhood
                        lastSolution.setInactive(v);
                    }
                }
                model.addConstr(fortSum >= 1);
                totalFortSize += forts[processedForts].size();
                if (output > 1) {
                    fmt::print("fort {:4}: {} #{}({})\n", processedForts, forts[processedForts], forts[processedForts].size(), blank);
                }
            }
            if (output) {
                fmt::print("#forts: {}, avg size: {:.2f}\n", forts.size(), double(totalFortSize) / double(forts.size()));
            }
            currentTime = now();
            remainingTimeout = std::max(1.0, timeLimit - std::chrono::duration_cast<std::chrono::seconds>(currentTime - startingTime).count());
            if (greedyUpper) {
                fastGreedy(lastSolution, (greedyUpper - 1) * 2);
                size_t initial = lastSolution.numActive();
                topDownGreedy(lastSolution, false, blankVertices);
                if (lastSolution.numActive() < upperBound) {
                    upperBound = lastSolution.numActive();
                    if (output) {
                        fmt::print("greedy {} → {}\n", initial, lastSolution.numActive());
                    }
                    std::swap(feasibleSolution, lastSolution);
                }
            }
            model.set(GRB_DoubleParam_TimeLimit, remainingTimeout);
            model.optimize();
            status = model.get(GRB_IntAttr_Status);
            size_t new_bound = 0;
            if (remainingTimeout <= 1.0) status = GRB_TIME_LIMIT;
            if (lastSolution.allObserved()) status = GRB_OPTIMAL;
            switch (status) {
                case GRB_INTERRUPTED:
                    if (static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound)) < upperBound) {
                        new_bound = std::max(lowerBound, static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound)));
                    }
                    break;
                case GRB_USER_OBJ_LIMIT:
                    new_bound = static_cast<size_t>(model.get(GRB_DoubleAttr_ObjVal));
                    break;
                case GRB_OPTIMAL:
                case GRB_TIME_LIMIT:
                    new_bound = static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound));
                    break;
                default:
                    fmt::print(stderr, "unexpected status: {}\n", status);
                    break;
            }
            for (auto v: state.graph().vertices()) {
                if (pi.at(v).get(GRB_DoubleAttr_X) > 0.5) {
                    lastSolution.setActive(v);
                } else if (state.isBlank(v)) {
                    lastSolution.setBlank(v);
                } else {
                    lastSolution.setInactive(v);
                }
            }
            if (lowerBound > new_bound) {
                if (status != GRB_TIME_LIMIT) { fmt::print("!!!lowerBound decreased {} → {}!!!\n", lowerBound, new_bound); }
            } else {
                lowerBound = new_bound;
            }
            if (lowerBound == upperBound) {
                std::swap(lastSolution, feasibleSolution);
                if (output) {
                    fmt::print("greedy solution is optimal\n");
                }
            }
            model.set(GRB_DoubleParam_BestObjStop, lowerBound);
            if(output) {
                fmt::print("LB: {}, UB: {} (status {})\n", lowerBound, upperBound, status);
            }
            if (status == GRB_USER_OBJ_LIMIT || status == GRB_INTERRUPTED) status = GRB_OPTIMAL;
            if (status != GRB_OPTIMAL || lastSolution.allObserved()) {
                callback(callback::When::FINAL, state, forts, lowerBound, upperBound);
                break;
            } else {
                callback(callback::When::INTERMEDIATE_HS, state, forts, lowerBound, upperBound);
            }
        }
        std::swap(state, lastSolution);
        if (greedyUpper) {
            fastGreedy(state, (greedyUpper - 1) * 2);
            upperBound = std::min(upperBound, state.numActive());
        } else {
            upperBound = state.numActive() + state.numBlank();
        }
        if (upperBound == lowerBound) {
            status = GRB_OPTIMAL;
        }

        switch (status) {
            case GRB_OPTIMAL:
                return {lowerBound, lowerBound, SolveState::Optimal};
            case GRB_TIME_LIMIT:
                return {lowerBound, upperBound, SolveState::Timeout};
            default:
                return {size_t{0}, state.numActive() + state.numBlank(), SolveState::Timeout};
        }
    } catch (const GRBException& ex) {
        fmt::print(stderr, "Gurobi Error [{}]: {}\n", ex.getErrorCode(), ex.getMessage());
        throw ex;
    }
}
} //namespace pds
