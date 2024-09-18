#include "fort_solve.hpp"

#include <gurobi_c++.h>
#include <range/v3/all.hpp>
#include <utility>

#include "pdssolve.hpp"
#include "gurobi_common.hpp"

namespace pds {

namespace {
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

struct Unimplemented {};
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

auto now() {
    return std::chrono::high_resolution_clock::now();
}

VertexList smithFortNeighborhood(const PdsState& state, bool output, double timeLimit, VertexSet& seen) {
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

std::vector<VertexList> initialFortsSmith(const PdsState& state, VertexSet& seen) {
    VertexSet junctions;
    std::vector<PowerGrid::VertexDescriptor> junctionVertices;
    for (auto v: state.graph().vertices()) {
        if (state.graph().degree(v) + state.isZeroInjection(v) >= 3) {
            junctions.insert(v);
            seen.insert(v);
            junctionVertices.push_back(v);
        }
    }
    auto paths = components(state, seen);
    for (auto v: junctionVertices) {
        seen.erase(v);
    }
    VertexMap<size_t> pathIndex;
    for (size_t i = 0 ; i < paths.size(); ++i) {
        for (auto v: paths[i]) {
            pathIndex.emplace(v, i);
        }
    }
    std::vector<VertexList> forts;
    for (auto v: state.graph().vertices()) {
        if (!pathIndex.contains(v)) {
            VertexList pathNeighbors;
            size_t type2path = -1;
            for (auto w: state.graph().neighbors(v)) {
                if (pathIndex.contains(w)) {
                    auto pathRep = paths[pathIndex[w]][0];
                    pathNeighbors.push_back(pathRep);
                    if (!seen.contains(pathRep)) {
                        seen.insert(pathRep);
                    } else {
                        type2path = pathRep;
                    }
                }
            }
            if (type2path < paths.size()) {
                VertexList fort;
                for (auto w: paths[type2path]) {
                    fort.push_back(w);
                }
                fort.push_back(v);
                forts.emplace_back(std::move(fort));
            }
            bool fortFound = false;
            std::vector<size_t> degreeOnePaths;
            for (auto rep: pathNeighbors) {
                for (auto u: paths[pathIndex[rep]]) {
                    bool hasOtherNeighbor = false;
                    for (auto w: state.graph().neighbors(u)) {
                        if (w != v && !pathIndex.contains(w)) {
                            hasOtherNeighbor = true;
                            if (!fortFound) {
                                for (auto z: state.graph().neighbors(w)) {
                                    if (pathIndex.contains(z)) {
                                        auto otherRep = paths[pathIndex[z]][0];
                                        if (otherRep != rep && seen.contains(otherRep)) {
                                            VertexList fort;
                                            for (auto t: paths[pathIndex[z]]) {
                                                fort.push_back(t);
                                            }
                                            for (auto t: paths[pathIndex[rep]]) {
                                                fort.push_back(t);
                                            }
                                            fort.push_back(v);
                                            forts.emplace_back(std::move(fort));
                                            fortFound = true;
                                        }
                                    }
                                }
                            } if (fortFound) break;
                        }
                    }
                    if (!hasOtherNeighbor) {
                        degreeOnePaths.push_back(pathIndex[rep]);
                    }
                    if (fortFound && degreeOnePaths.size() > 1) break;
                } if (fortFound && degreeOnePaths.size() > 1) break;
            }
            if (degreeOnePaths.size() > 1) {
                VertexList fort;
                for (auto p: degreeOnePaths) {
                    for (auto w: paths[p]) {
                        fort.push_back(w);
                    }
                }
                fort.push_back(v);
                forts.emplace_back(std::move(fort));
            }
            for (auto w: pathNeighbors) {
                seen.erase(w);
            }
        }
    }

    fmt::print("found {} initial forts\n", forts.size());
    return forts;
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
    auto inactive = state.graph().vertices()
                 | ranges::views::filter([&state](auto v) { return state.isInactive(v); })
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
    for (auto v: inactive) { state.setInactive(v); }
    return forts;
}


std::vector<VertexList> initializeForts(PdsState& state, int variant, VertexSet& seen) {
    switch (variant) {
        case 1:
            return initialForts(state, seen);
        case 2:
            return initialForts2(state, seen);
        case 3:
            return initialForts3(state, seen);
        case 4:
            return initialFortsSmith(state, seen);
        case 0:
        default:
            return {};
    }
}
auto violatedForts(PdsState& lastSolution, int variant, double remainingTimeout, VertexSet& seen, int output) -> std::vector<VertexList> {
    switch(variant) {
        case 1:
            return {smithFortNeighborhood(lastSolution, false, remainingTimeout, seen)};
        case 2:
            return {bozemanFortNeighborhood2(lastSolution, false, remainingTimeout)};
        case 3:
            return {bozemanFortNeighborhood3(lastSolution, false, remainingTimeout)};
        case 4: {
            if (output) { fmt::print("finding forts in {} blank vertices\n", lastSolution.numBlank()); }
            return initialForts3(lastSolution, seen);
        }
        case 0:
        default:
            return {bozemanFortNeighborhood(lastSolution, output, remainingTimeout, seen)};
    }
}
struct Callback : public GRBCallback {
    size_t* lower;
    size_t* upper;
    VertexMap<GRBVar>* pi;
    const PdsState* base;
    PdsState* solution;
    PdsState* upperBound;
    std::span<PdsState::Vertex> blank;
    int earlyStop;
    int intermediateForts;
    std::vector<VertexList>* forts;
    VertexSet* seen;
    int output;
    Callback(size_t& lower, size_t& upper, VertexMap<GRBVar>& pi,
             const PdsState& base, PdsState& solution, PdsState & upperBound,
             std::span<PdsState::Vertex> blank, int earlyStop,
             int intermediateForts, std::vector<VertexList>& forts, VertexSet& seen,
             int output)
            : lower(&lower), upper(&upper), pi(&pi), base(&base), solution(&solution),
              upperBound(&upperBound), blank(blank), earlyStop(earlyStop),
              intermediateForts(intermediateForts), forts(&forts), seen(&seen), output(output)
    { }
    void callback() override {
        if (where == GRB_CB_MIP) {
            if (getIntInfo(GRB_CB_MIP_SOLCNT) > 0) {
                auto objVal = static_cast<size_t>(getDoubleInfo(GRB_CB_MIP_OBJBST) + 0.5);
                auto objBound = static_cast<size_t>(getDoubleInfo(GRB_CB_MIP_OBJBND));
                if (objVal <= *upper && !solution->allObserved() && earlyStop > 1 && (objBound > *lower || earlyStop > 2)) {
                    abort();
                }
            }
        }
        if (where == GRB_CB_MIPSOL) {
            auto objVal = static_cast<size_t>(getDoubleInfo(GRB_CB_MIPSOL_OBJ) + 0.5);
            auto objBound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
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

                if (!solution->allObserved()) {
                    if (intermediateForts >= 0) {
                        for (auto& f: violatedForts(*solution, intermediateForts, 10, *seen, false)) {
                            forts->push_back(std::move(f));
                        }
                    }

                    if (output > 1) {
                        fmt::print("new solution {} <= {}; {} <= {}; unobserved #{}\n", *lower, objBound, objVal, *upper, solution->numUnobserved());
                    }

                    if (earlyStop > 1 && objBound > *lower) { abort(); }
                } else if (objVal < *upper) {
                    *upper = objVal;
                    *upperBound = *solution;
                    fmt::print("new power dominating set {} <= {}; {} <= {}\n", *lower, objBound, objVal, *upper, objVal);
                    //if (*upper <= *lower && earlyStop) { abort(); }
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
        int fortInit,
        int greedyUpper,
        int earlyStop,
        callback::FortCallback fortCallback,
        BoundCallback boundCallback,
        int intermediateForts
) {
    auto lastSolution = state;
    //writePds(lastSolution.graph(), fmt::format("out/0_input.pds"));
    auto blankVertices = state.graph().vertices()
                         | ranges::views::filter([&state] (auto v) { return state.isBlank(v); })
                         | ranges::to<std::vector<PdsState::Vertex>>;
    PdsState feasibleSolution = state;

    size_t lowerBound = 0;
    size_t upperBound = state.numActive() + state.numBlank();
    if (state.allObserved()) {
        return {state.numActive(), state.numActive(), SolveState::Optimal};
    } else if (state.numBlank() == 0) {
        return {state.graph().numVertices(), 0, SolveState::Infeasible};
    }
    try {
        VertexSet seen;
        auto forts = initializeForts(state, fortInit, seen);
        GRBModel model(getEnv());
        model.set(GRB_IntParam_LogToConsole, false);
        model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        model.set(GRB_IntParam_NumericFocus, 1);
        model.set(GRB_DoubleParam_MIPGap, 1e-6);
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
        Callback cb(
                lowerBound, upperBound,
                pi, state, lastSolution, feasibleSolution,
                blankVertices, earlyStop, intermediateForts, forts, seen, output);
        model.setCallback(&cb);
        auto startingTime = now();
        int status = 0;
        size_t processedForts = 0;
        size_t totalFortSize = 0;
        while (true) {
            auto currentTime = now();
            double remainingTimeout = std::max(0.0, timeLimit - std::chrono::duration_cast<std::chrono::duration<double>>(currentTime - startingTime).count());
            auto moreForts = violatedForts(lastSolution, variant, remainingTimeout, seen, output);
            for (auto f: moreForts) {
                forts.emplace_back(std::move(f));
            }
            if (forts.empty()) {
                return {state.graph().numVertices(), 0, SolveState::Infeasible};
            }
            if (forts.back().empty()) break;
            for (; processedForts < forts.size(); ++processedForts) {
                GRBLinExpr fortSum;
                size_t blank = 0;
                for (auto& v: forts[processedForts]) {
                    if (!state.graph().hasVertex(v)) {
                        fmt::print("!!!invalid vertex {} in fort: {}!!!\n", v, forts[processedForts]);
                    }
                    if (state.isActive(v)) {
                        if (output) {
                            fmt::print("???active vertex in neighborhood??? {} {} {}\n",
                                       v,
                                       lastSolution.numObserved(),
                                       lastSolution.graph().numVertices());
                        }
                    }
                    if (state.isBlank(v)) {
                        fortSum += pi.at(v);
                        ++blank;
                    }
                }
                if (blank == 0) {
                    fmt::print("??? infeasible fort: {} has no blank vertices\n", forts[processedForts]);
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
            remainingTimeout = std::max(0.0, timeLimit - std::chrono::duration_cast<std::chrono::duration<double>>(currentTime - startingTime).count());
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
            switch (status) {
                case GRB_INFEASIBLE:
                    return {1, 0, SolveState::Infeasible};
                case GRB_INTERRUPTED:
                    if (static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound)) < state.graph().numVertices()) {
                        new_bound = std::max(lowerBound, static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound)));
                    }
                    break;
                case GRB_USER_OBJ_LIMIT:
                case GRB_OPTIMAL:
                    new_bound = static_cast<size_t>(model.get(GRB_DoubleAttr_ObjVal) + 0.05);
                    break;
                case GRB_TIME_LIMIT:
                    new_bound = static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound));
                    break;
                default:
                    fmt::print(stderr, "unexpected status: {}\n", status);
                    break;
            }
            if (model.get(GRB_IntAttr_SolCount) > 0) {
                for (auto v: state.graph().vertices()) {
                    if (pi.at(v).get(GRB_DoubleAttr_X) > 0.5) {
                        lastSolution.setActive(v);
                    } else if (state.isBlank(v)) {
                        lastSolution.setBlank(v);
                    } else {
                        lastSolution.setInactive(v);
                    }
                }
            }
            if (lowerBound > new_bound) {
                if (status != GRB_TIME_LIMIT) {
                    fmt::print("!!!lowerBound decreased {} → {}!!!\n", lowerBound, new_bound);
                    if (size_t(model.get(GRB_DoubleAttr_ObjVal))) {
                        fmt::print("!!!wrong lower bound: obj {} < bound {}!!!\n", size_t(GRB_DoubleAttr_ObjVal), lowerBound);
                    }
                }
            } else if (new_bound <= state.graph().numVertices()) { // only use valid bounds
                lowerBound = new_bound;
            }
            if (lowerBound == upperBound) {
                std::swap(lastSolution, feasibleSolution);
                if (output) {
                    fmt::print("greedy solution is optimal\n");
                }
            } else if (lastSolution.allObserved()) {
                upperBound = lastSolution.numActive();
            }
            boundCallback(lowerBound, upperBound, forts.size());
            if(earlyStop > 0) {
                model.set(GRB_DoubleParam_BestObjStop, double(lowerBound));
            }
            if(output) {
                fmt::print("LB: {}, UB: {} (status {}) (local LB: {}, UP: {})\n", lowerBound, upperBound, status, model.get(GRB_DoubleAttr_ObjBound), model.get(GRB_DoubleAttr_ObjVal));
            }
            if (remainingTimeout <= 1.0) status = GRB_TIME_LIMIT;
            if (lastSolution.allObserved()) {
                feasibleSolution = lastSolution;
                fortCallback(callback::When::FINAL, state, forts, lowerBound, upperBound);
                break;
            } else {
                fortCallback(callback::When::INTERMEDIATE_HS, state, forts, lowerBound, upperBound);
            }
            for (auto v: feasibleSolution.graph().vertices()) {
                if (feasibleSolution.isActive(v)) {
                    pi.at(v).set(GRB_DoubleAttr_Start, 1.0);
                } else {
                    pi.at(v).set(GRB_DoubleAttr_Start, 0.0);
                }
            }
        }
        std::swap(state, feasibleSolution);
        if (upperBound == lowerBound) {
            status = GRB_OPTIMAL;
        }
        if (upperBound < state.numActive()) {
            fmt::print("!!!upper bound too low!!!\n");
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

void addFortConstraints(MIPModel &mipmodel, PdsState &state, int fortInit) {
    auto& model = *mipmodel.model;
    const auto& xi = mipmodel.xi;
    VertexSet seen;
    auto forts = initializeForts(state, fortInit, seen);
    for (auto & fort: forts) {
        GRBLinExpr fortSum;
        for (auto v: fort) {
            if (xi.contains(v)) {
                fortSum += xi.at(v);
            } else {
                fmt::print("!!! vertex not in mip model: {}!!!\n", v);
            }
        }
        model.addConstr(fortSum >= 1);
    }
}
namespace {
struct LazyFortCB : public GRBCallback {
    PdsState* solution;
    PdsState bestSolution;
    int output;
    callback::FortCallback fortCB;
    BoundCallback boundCB;
    VertexList blank;
    VertexSet initialBlank;
    std::vector<VertexList> forts;
    VertexMap<GRBVar> xi;
    VertexSet seen;
    GRBModel model;
    size_t lower;
    LazyFortCB(PdsState& input, int output, double timeLimit, callback::FortCallback fortCB, BoundCallback boundCB) :
            solution(&input), bestSolution(input), output(output), fortCB(std::move(fortCB)), boundCB(std::move(boundCB)),
            seen{}, model(getEnv()), lower(0)
    {
        model.set(GRB_IntParam_LogToConsole, output >= 2);
        model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        model.set(GRB_IntParam_LazyConstraints, 1);
        //model.set(GRB_IntParam_NumericFocus, 1);
        //model.set(GRB_DoubleParam_MIPGap, 1e-6);
        model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
        model.setCallback(this);
        fastGreedy(bestSolution, false);
        for (auto v: input.graph().vertices()) {
            xi.emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("x_{}", v)));
            if (input.isBlank(v)) {
                blank.push_back(v);
                initialBlank.insert(v);
            } else if (input.isInactive(v)) {
                model.addConstr(xi.at(v) == 0);
            } else if (input.isActive(v)) {
                model.addConstr(xi.at(v) == 1);
            }
        }
        auto newForts = initialForts3(*solution, seen);
        for (auto &f: newForts) {
            GRBLinExpr fortSum;
            for (auto v: f) {
                fortSum += xi.at(v);
            }
            model.addConstr(fortSum >= 1);
            forts.emplace_back(std::move(f));
        }
        //newForts = initialForts3(*solution, seen);
        //for (auto &f: newForts) {
        //    GRBLinExpr fortSum;
        //    for (auto v: f) {
        //        fortSum += xi.at(v);
        //    }
        //    model.addConstr(fortSum >= 1);
        //    forts.emplace_back(std::move(f));
        //}
        lower = 0;
        this->boundCB(lower, bestSolution.numActive(), 0);
    }
    SolveResult solve() {
        model.optimize();
        std::swap(*solution, bestSolution);
        switch (model.get(GRB_IntAttr_Status)) {
            case GRB_OPTIMAL:
                return SolveResult{solution->numActive(), solution->numActive(), SolveState::Optimal};
            case GRB_TIME_LIMIT:
                return SolveResult{static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound)), solution->numActive(), SolveState::Timeout};
            case GRB_INFEASIBLE:
                return SolveResult{solution->graph().numVertices(), 0, SolveState::Infeasible};
            default:
                return SolveResult{0, solution->graph().numVertices(), SolveState::Other};
        }
    }

    void callback() override {
        switch (where) {
            case GRB_CB_MIPSOL: {
                for (auto v: blank) {
                    if (getSolution(xi.at(v)) > 0.5) {
                        solution->setActive(v);
                    } else if (initialBlank.contains(v)) {
                        solution->setBlank(v);
                    } else {
                        fmt::print("!!!changed inactive vertex {}!!!\n", v);
                        solution->setInactive(v);
                    }
                }
                if (!solution->allObserved()) {
                    addLazyForts();
                    addLazyForts();
                    addLazyForts();
                    lower = std::max(lower, getLower(getDoubleInfo(GRB_CB_MIPSOL_OBJBND)));
                    // compute new upper actual bound
                    fastGreedy(*solution, false);
                    topDownGreedy(*solution, false, blank);
                    fortCB(callback::When::INTERMEDIATE_HS, *solution, forts, lower, bestSolution.numActive());
                }
                if (solution->numActive() < bestSolution.numActive()) {
                    bestSolution = *solution;
                    boundCB(lower, bestSolution.numActive(), forts.size());
                }
                if (output > 0) {
                    fmt::print("LB: {}\tUB: {}\t#F: {}\n", lower, bestSolution.numActive(), forts.size());
                }
                break;
            }
            case GRB_CB_MIP: {
                size_t newLower = getLower(getDoubleInfo(GRB_CB_MIP_OBJBND));
                if (newLower > lower) {
                    lower = newLower;
                    boundCB(lower, bestSolution.numActive(), forts.size());
                    if (output > 0) {
                        fmt::print("LB: {}\tUB: {}\t#F: {}\n", lower, bestSolution.numActive(), forts.size());
                    }
                }
                break;
            }
        }
    }
private:
    size_t getLower(double grbLower) {
        if (grbLower >= GRB_INFINITY) {
            return 0;
        } else {
            return static_cast<size_t>(std::max(0.0, grbLower));
        }
    }
    void addLazyForts() {
        auto newForts = initialForts3(*solution, seen);
        for (auto &f: newForts) {
            GRBLinExpr fortSum;
            for (auto v: f) {
                fortSum += xi.at(v);
            }
            addLazy(fortSum >= 1);
            forts.emplace_back(std::move(f));
        }
    }
};
}

SolveResult solveLazyForts(PdsState &state, int output, double timeLimit, callback::FortCallback fortCB, BoundCallback boundsCB) {
    LazyFortCB lazyForts(state, output, timeLimit, std::move(fortCB), std::move(boundsCB));
    return lazyForts.solve();
}
} //namespace pds
