//
// Created by max on 10.08.22.
//
#include "gurobi_solve.hpp"

#include <gurobi_c++.h>
#include <unistd.h>

#include "graphio.hpp"

namespace pds {

namespace {
GRBEnv &getEnvWithOutput() {
    static thread_local GRBEnv env;
    return env;
}
GRBEnv &getEnv() {
    // Suppress Academic License Info by redicreting stdout
    static thread_local bool loaded = false;
    if (!loaded) {
        auto oldstdout = dup(STDOUT_FILENO);
        FILE *result = freopen("/dev/null", "a", stdout);
        if (!result) throw std::ios_base::failure(strerror(errno), std::error_code(errno, std::system_category()));
        auto& env = getEnvWithOutput();
        dup2(oldstdout, STDOUT_FILENO);
        loaded = true;
        return env;
    } else {
        return getEnvWithOutput();
    }
}
}

MIPModel::MIPModel() : model(std::make_unique<GRBModel>(getEnv())), xi() { }
MIPModel::~MIPModel() { }

void preloadMIPSolver() {
    getEnv();
}

SolveState solveModel(PdsState& state, map<PowerGrid::VertexDescriptor, GRBVar>& xi, GRBModel& model) {
    model.optimize();
    if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
        return SolveState::Infeasible;
    } else {
        //std::cout << "result: " << model.get(GRB_DoubleAttr_ObjBound) << std::endl;
        for (auto v: state.graph().vertices()) {
            if (xi.at(v).get(GRB_DoubleAttr_X) > 0.5) {
                state.setActive(v);
            } else {
                state.setInactive(v);
            }
        }
    }
    switch (model.get(GRB_IntAttr_Status)) {
    case GRB_OPTIMAL:
        return SolveState::Optimal;
    case GRB_TIME_LIMIT:
        return SolveState::Timeout;
    default:
        return SolveState::Other;
    }
}

SolveResult solveMIP(MIPModel & mipmodel, bool output, double timeLimit) {
    auto& model = *mipmodel.model;
    auto oldstdout = dup(STDOUT_FILENO);
    FILE *result = freopen("/dev/null", "a", stdout);
    if (!result) throw std::ios_base::failure(strerror(errno), std::error_code(errno, std::system_category()));
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    dup2(oldstdout, STDOUT_FILENO);
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);

    model.optimize();
    switch (model.get(GRB_IntAttr_Status)) {
        case GRB_INFEASIBLE:
            return { size_t{1}, size_t{0}, SolveState::Infeasible };
        case GRB_OPTIMAL:
            return { size_t(model.get(GRB_DoubleAttr_ObjBound)), size_t(model.get(GRB_DoubleAttr_ObjVal)), SolveState::Optimal };
        case GRB_TIME_LIMIT:
            return { size_t(model.get(GRB_DoubleAttr_ObjBound)), size_t(model.get(GRB_DoubleAttr_ObjVal)), SolveState::Timeout };
        default:
            return { size_t{1}, size_t{0}, SolveState::Other };
    }
}

void applySolution(PdsState& state, MIPModel &model) {
    auto& xi = model.xi;
    for (auto v: state.graph().vertices()) {
        if (xi.at(v).get(GRB_DoubleAttr_X) > 0.5) {
            state.setActive(v);
        } else {
            state.setInactive(v);
        }
    }
}

MIPModel modelJovanovicExpanded(PdsState& state) {
    namespace r3 = ranges;
    try {
        MIPModel mipmodel;
        auto& model = *mipmodel.model;
        auto& xi = mipmodel.xi;
        map<PowerGrid::VertexDescriptor, GRBVar> si;
        map<std::pair<PowerGrid::VertexDescriptor, PowerGrid::VertexDescriptor>, GRBVar> pij;
        auto &graph = static_cast<const PdsState &>(state).graph();
        const double M = 2 * graph.numVertices();
        {
            for (auto v: graph.vertices()) {
                xi.insert({v, model.addVar(0.0, M, 1.0, GRB_BINARY)});
                // si[v] >= 1 && si[v] <= num_vertices(graph)
                si.insert({v, model.addVar(1.0, M, 0.0, GRB_CONTINUOUS)});
            }
            for (auto e: graph.edges()) {
                auto v = graph.source(e);
                auto w = graph.target(e);
                if (graph[v].zero_injection)
                    pij[{v, w}] = model.addVar(0.0, M, 0.0, GRB_BINARY); //pij_p[2 * i];
                if (graph[w].zero_injection)
                    pij[{w, v}] = model.addVar(0.0, M, 0.0, GRB_BINARY); //pij_p[2 * i + 1];
            }
            //model.setObjective(objective, GRB_MINIMIZE);
        }

        for (auto v: graph.vertices()) {
            model.addConstr(si.at(v) >= 1);
            for (auto w: r3::concat_view(graph.neighbors(v), r3::single_view{v})) {
                switch (state.activeState(w)) {
                    case PmuState::Blank:
                        model.addConstr(si.at(v) <= xi.at(w) + M * (1 - xi.at(w)));
                        break;
                    case PmuState::Inactive:
                        model.addConstr(si.at(v) <= M);
                        break;
                    case PmuState::Active:
                        model.addConstr(si.at(v) <= 1);
                }
            }
            GRBLinExpr observingNeighbors = xi.at(v);
            // v observes at most one neighbor w1
            GRBLinExpr outObserved;
            // v is observed by at most one neighbor w1
            GRBLinExpr inObserved;
            for (auto w: graph.neighbors(v)) {
                observingNeighbors += xi[w];
                if (graph[w].zero_injection) {
                    observingNeighbors += pij[{w, v}];
                    inObserved += pij[{w, v}];
                }
                if (graph[v].zero_injection) {
                    outObserved += pij[{v, w}];
                }
            }
            if (!state.isActive(v)) {
                model.addConstr(si[v] <= M * observingNeighbors);
            }
            model.addConstr(inObserved <= !state.isActive(v));
            if (graph[v].zero_injection) {
                model.addConstr(outObserved <= 1);
            }

            model.addConstr(si[v] <= graph.numVertices());

            if (state.isActive(v)) {
                model.addConstr(xi[v] == 1);
            } else if (state.isInactive(v)) {
                model.addConstr(xi[v] == 0);
            }
        }

        for (auto e: graph.edges()) {
            auto v = graph.source(e);
            auto w = graph.target(e);
            //model.addConstr(pij[{v, w}] + pij[{w, v}] <= 1);
            for (auto t: r3::concat_view(graph.neighbors(w), r3::single_view{w})) {
                if (t != v) {
                    if (graph[w].zero_injection) {
                        model.addConstr(si[v] >= si[t] + 1 - M * (1 - pij[{w, v}]));
                    }
                }
            }
            for (auto t: r3::concat_view(graph.neighbors(v), r3::single_view{v})) {
                if (t != w) {
                    if (graph[v].zero_injection) {
                        model.addConstr(si[w] >= si[t] + 1 - M * (1 - pij[{v, w}]));
                    }
                }
            }
        }

        return mipmodel;
    } catch (GRBException ex) {
        fmt::print(stderr, "Gurobi Exception {}: {}\n", ex.getErrorCode(), ex.getMessage());
        throw ex;
    }
}

MIPModel modelJovanovic(PdsState& state) {
    MIPModel mipmodel;
    auto& model = *mipmodel.model;
    auto& xi = mipmodel.xi;
    map<PowerGrid::VertexDescriptor, GRBVar> si;
    map<std::pair<PowerGrid::VertexDescriptor, PowerGrid::VertexDescriptor>, GRBVar> pij;
    auto &graph = static_cast<const PdsState &>(state).graph();
    const auto N = static_cast<double>(graph.numVertices());
    const auto M = static_cast<double>(2 * graph.numVertices());
    for (auto v: graph.vertices()) {
        xi.try_emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY));
        xi[v].set(GRB_StringAttr_VarName, fmt::format("x_{}", v));
        si.try_emplace(v, model.addVar(1.0, N + 1, 0.0, GRB_CONTINUOUS));
        si[v].set(GRB_StringAttr_VarName, fmt::format("s_{}", v));
        for (auto w: graph.neighbors(v)) {
            pij.try_emplace({v, w}, model.addVar(0.0, 1.0, 0.0, GRB_BINARY));
            pij[{v, w}].set(GRB_StringAttr_VarName, fmt::format("p_{},{}", v, w));
            if (!state.isZeroInjection(v)) {
                model.addConstr(pij[{v, w}] == 0.0);
            }
        }
        if (state.isActive(v)) {
            model.addConstr(xi[v] == 1.0);
        }
        if (state.isInactive(v)) {
            model.addConstr(xi[v] == 0.0);
        }
    }
    for (auto v: graph.vertices()) {
        //model.addConstr(si[v] <= xi[v] + M * (1.0 - xi[v]));
        GRBLinExpr outSum;
        GRBLinExpr inSum;
        GRBLinExpr observingNeighbors = xi.at(v);
        GRBLinExpr observingEdges;
        for (auto w: graph.neighbors(v)) {
            //model.addConstr(si[v] <= xi[w] + M * (1.0 - xi[w]));
            outSum += pij[{v, w}];
            inSum += pij[{w, v}];
            observingNeighbors += xi[w];
            observingEdges += pij[{w, v}];
        }
        model.addConstr(si[v] <= M * (observingNeighbors + observingEdges));
        model.addConstr(outSum <= 1);
        model.addConstr(inSum <= 1);
        model.addConstr(si[v] <= N);
        for (auto w: graph.neighbors(v)) {
            for (auto t: ranges::concat_view(graph.neighbors(w), ranges::single_view(w))) {
                if (t != v) {
                    model.addConstr(si[v] >= si[t] + 1 - M * (1.0 - pij[{w, v}]));
                }
            }
        }
    }
    return mipmodel;
}

MIPModel modelDomination(PdsState& state) {
    MIPModel mipmodel;
    auto& model = *mipmodel.model;
    auto& xi = mipmodel.xi;
    for (auto v: state.graph().vertices()) {
        xi.try_emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY));
    }
    for (auto v: state.graph().vertices()) {
        GRBLinExpr sum = xi.at(v);
        for (auto w: state.graph().neighbors(v)) {
            sum += xi.at(w);
        }
        model.addConstr(sum >= 1);
    }

    return mipmodel;
}

MIPModel modelBrimkovExpanded(PdsState& state) {
    MIPModel mipmodel;
    auto& model = *mipmodel.model;
    auto& xi = mipmodel.xi;
    map <PowerGrid::VertexDescriptor, GRBVar> si;
    map <PowerGrid::EdgeDescriptor, GRBVar> ye;
    auto T = static_cast<double>(state.graph().numVertices());
    for (auto v: state.graph().vertices()) {
        xi.try_emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, "x"));
        si.try_emplace(v, model.addVar(0.0, static_cast<double>(T), 0.0, GRB_INTEGER, "s"));
        for (auto e: state.graph().outEdges(v)) {
            assert(!ye.contains(e));
            ye.try_emplace(e, model.addVar(0.0, 1.0, 0.0, GRB_BINARY));
        }
    }
    for (auto u: state.graph().vertices()) {
        // x_v + \sum_{e \in N(v)} y_e = 1 (3)
        if (state.isActive(u)) {
            model.addConstr(xi.at(u) == 1);
        } else {
            GRBLinExpr observers = 0;
            if (state.isBlank(u)) observers += xi.at(u);
            for (auto e: state.graph().inEdges(u)) {
                auto v = state.graph().source(e);
                assert(v != u);
                if (!state.isInactive(v)) observers += xi.at(v);
                if (state.isZeroInjection(v)) observers += ye.at(e);
            }
            model.addConstr(observers >= 1);
        }
        if (!state.isZeroInjection(u) && !state.isInactive(u)) {
            GRBLinExpr outPropagation;
            for (auto e: state.graph().outEdges(u)) {
                outPropagation += ye.at(e);
            }
            model.addConstr(outPropagation <= 0);
        } else {
            GRBLinExpr outPropagation;
            for (auto e: state.graph().outEdges(u)) {
                outPropagation += ye.at(e);
            }
            model.addConstr(outPropagation <= 1);
        }
        for (auto e: state.graph().outEdges(u)) {
            auto v = state.graph().target(e);
            {
                // s_u - s_v + (T + 1) y_{uv} <= T (4) e@(u,v) \in E'
                model.addConstr(si.at(u) - si.at(v) + (T + 1) * ye.at(e) <= T);
                for (auto w: state.graph().neighbors(u)) {
                    if (v != w) {
                        unused(u, w, v, e);
                        // s_w - s_v + (T + 1) y_e <= T + (T+1) x_u , e@(u,v) \in E', w \in N(u) - v
                        if (state.isBlank(u)) {
                            if (state.isZeroInjection(u)) {
                                model.addConstr(si.at(w) - si.at(v) + (T + 1) * ye.at(e) <= T + (T + 1) * xi.at(u));
                            }
                        } else if (state.isActive(u)) {
                            // ignore, no action needed
                        } else if (state.isInactive(u)) {
                            if (state.isZeroInjection(u)) {
                                model.addConstr(si.at(w) - si.at(v) + (T + 1) * ye.at(e) <= T);
                            }
                        }
                    }
                }
            }
        }
    }
    return mipmodel;
}

MIPModel modelBrimkov(PdsState& state) {
    MIPModel mipmodel;
    auto& model = *mipmodel.model;
    auto& xi = mipmodel.xi;
    map <PowerGrid::VertexDescriptor, GRBVar> si;
    map <PowerGrid::EdgeDescriptor, GRBVar> ye;
    size_t T = state.graph().numVertices();
    for (auto v: state.graph().vertices()) {
        xi.try_emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, "x"));
        si.try_emplace(v, model.addVar(0.0, static_cast<double>(T), 0.0, GRB_INTEGER, "s"));
        for (auto e: state.graph().outEdges(v)) {
            assert(!ye.contains(e));
            ye.try_emplace(e, model.addVar(0.0, 1.0, 0.0, GRB_BINARY));
        }
    }
    for (auto u: state.graph().vertices()) {
        // x_v + \sum_{e \in N(v)} y_e = 1 (3)
        GRBLinExpr observers = 0;
        observers += xi.at(u);
        for (auto e: state.graph().inEdges(u)) {
            auto v = state.graph().source(e);
            assert(v != u);
            observers += ye.at({v, u});
        }
        model.addConstr(observers == 1);
        for (auto e: state.graph().outEdges(u)) {
            auto v = state.graph().target(e);
            if (!state.isZeroInjection(u)) {
                model.addConstr(si.at(v) <= si.at(u) + 1);
            }
         // s_u - s_v + (T + 1) y_{uv} <= T (4) e@(u,v) \in E'
            model.addConstr(si.at(u) - si.at(v) + (T + 1) * ye.at(e) <= T);
            for (auto w: state.graph().neighbors(u)) {
                if (v != w) {
                    unused(u, w, v, e);
         // s_w - s_v + (T + 1) y_e <= T + (T+1) x_u , e@(u,v) \in E', w \in N(u) - v
                    model.addConstr(si.at(w) - si.at(v) + (T + 1) * ye.at(e) <= T + (T+1) * xi.at(u));
                }
            }
        }
    }
    return mipmodel;
}

MIPModel modelAzamiBrimkov(PdsState& state) {
    MIPModel mipmodel;
    auto& model = *mipmodel.model;
    auto &graph = static_cast<const PdsState &>(state).graph();
    auto& xi = mipmodel.xi;
    map <std::pair<size_t, PowerGrid::VertexDescriptor>, GRBVar> si;
    map <std::pair<PowerGrid::VertexDescriptor, PowerGrid::VertexDescriptor>, GRBVar> ye;
    size_t T = graph.numVertices();
    for (auto v: graph.vertices()) {
        xi.try_emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY));
        GRBLinExpr vObserved;
        for (size_t i = 0; i < T; ++i) {
            si.try_emplace({i, v}, model.addVar(0.0, 1.0, 0.0, GRB_BINARY));
            vObserved += si[{i, v}];
        }
        model.addConstr(vObserved >= 0);
        for (auto w: graph.neighbors(v)) {
            ye.try_emplace({v, w}, model.addVar(0.0, 1.0, 0.0, GRB_BINARY));
        }
    }
    struct Unimplemented{ };
    throw Unimplemented{};
    return mipmodel;
}

using VertexList = std::vector<PowerGrid::VertexDescriptor>;

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

SolveResult solveBozeman(PdsState &state, int output, double timeLimit, int variant) {
    unused(state, output, timeLimit);
    auto lastSolution = state;
    //writePds(lastSolution.graph(), fmt::format("out/0_input.pds"));

    size_t bound = 0;
    try {
        VertexSet seen;
        auto forts = initialForts2(state, seen);
        GRBModel model(getEnv());
        model.set(GRB_IntParam_LogToConsole, false);
        model.set(GRB_DoubleParam_TimeLimit, timeLimit);
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
        auto startingTime = now();
        int status;
        size_t processedForts = 0;
        size_t solvedNeighborhoods = 0;
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
                    if (output && lastSolution.isActive(v)) {
                        fmt::print("!!!active vertex in neighborhood!!! {} {} {}\n",
                                   v,
                                   lastSolution.numObserved(),
                                   lastSolution.graph().numVertices());
                    }
                    if (state.isBlank(v)) {
                        fortSum += pi.at(v);
                        ++blank;
                    }
                }
                model.addConstr(fortSum >= 1);
                totalFortSize += forts[processedForts].size();
                if (output > 1) {
                    fmt::print("fort {:4}: {} #{}({})\n", processedForts, forts[processedForts], forts[processedForts].size(), blank);
                }
            }
            if (output) {
                fmt::print("#forts: {}, avg size: {:.2f}\n", forts.size(), double(totalFortSize) / forts.size());
            }
            currentTime = now();
            remainingTimeout = std::max(1.0, timeLimit - std::chrono::duration_cast<std::chrono::seconds>(currentTime - startingTime).count());
            {
                auto filename = fmt::format("hs/{:04}.hs", solvedNeighborhoods);
                ++solvedNeighborhoods;
                FILE* file = fopen(filename.c_str(), "w");
                fmt::print(file, "{} {}\n", state.graph().numVertices(), forts.size());
                for (auto& f: forts) {
                    VertexList fortBlank;
                    for (auto v: f) {
                        if (state.isBlank(v)) {
                            fortBlank.push_back(v);
                        }
                    }
                    ranges::sort(fortBlank);
                    fmt::print(file, "{}\n", fmt::join(fortBlank, " "));
                }
                fclose(file);
            }
            model.set(GRB_DoubleParam_TimeLimit, remainingTimeout);
            model.optimize();
            status = model.get(GRB_IntAttr_Status);
            size_t new_bound = 0;
            if (remainingTimeout <= 1.0) status = GRB_TIME_LIMIT;
            if (lastSolution.allObserved()) status = GRB_OPTIMAL;
            switch (status) {
                case GRB_USER_OBJ_LIMIT:
                    new_bound = model.get(GRB_DoubleAttr_ObjVal);
                    break;
                case GRB_OPTIMAL:
                case GRB_TIME_LIMIT:
                    new_bound = model.get(GRB_DoubleAttr_ObjBound);
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
            if (bound > new_bound) { fmt::print(stderr, "!!!bound decreased!!!\n"); }
            bound = new_bound;
            model.set(GRB_DoubleParam_BestObjStop, bound);
            if(output) {
                fmt::print("LB: {} (status {})\n", bound, status);
            }
            if (status == GRB_USER_OBJ_LIMIT) status = GRB_OPTIMAL;
            if (status != GRB_OPTIMAL || lastSolution.allObserved()) break;
        }
        std::swap(state, lastSolution);

        switch (status) {
            case GRB_OPTIMAL:
                return {bound, bound, SolveState::Optimal};
            case GRB_TIME_LIMIT:
                return {bound, (state.numActive() + state.numBlank()), SolveState::Timeout};
            default:
                return {size_t{0}, state.numActive() + state.numBlank(), SolveState::Timeout};
        }
    } catch (const GRBException& ex) {
        fmt::print(stderr, "Gurobi Error [{}]: {}\n", ex.getErrorCode(), ex.getMessage());
        throw ex;
    }

}

} // namespace pds
