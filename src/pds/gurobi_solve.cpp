//
// Created by max on 10.08.22.
//
#include "gurobi_solve.hpp"

#include <gurobi_c++.h>

namespace pds {

namespace {
GRBEnv &getEnv() {
    // Suppress Academic License Info by redicreting stdout
    auto oldstdout = dup(STDOUT_FILENO);
    FILE* result = freopen("/dev/null", "a", stdout);
    if (!result) throw std::ios_base::failure(strerror(errno), std::error_code(errno, std::system_category()));
    static thread_local GRBEnv env;
    dup2(oldstdout, STDOUT_FILENO);
    return env;
}
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

//SolveState solve_pds(PdsState& state, bool output, double timeLimit) {
//    auto active = state.active();
//    auto observed = state.observed();
//    bool result = solve_pds(state.graph(), active, output, timeLimit);
//    for (auto v: state.graph().vertices()) {
//        if (active[v] == PmuState::Active) {
//            state.setActive(v);
//        }
//    }
//    return result;
//}

SolveState solve_pds(PdsState& state, bool output, double timeLimit) {
    namespace r3 = ranges;
    try {
        auto model = GRBModel(getEnv());
        model.set(GRB_IntParam_LogToConsole, int{output});
        model.set(GRB_StringParam_LogFile, "gurobi.log");
        model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        //GRBVar *xi_p = model.addVars(num_vertices(graph), GRB_BINARY);
        //GRBVar *pij_p = model.addVars(2 * num_edges(graph), GRB_BINARY);
        //GRBVar *si_p = model.addVars(num_vertices(graph));
        map<PowerGrid::VertexDescriptor, GRBVar> xi;
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

        return solveModel(state, xi, model);
    } catch (GRBException ex) {
        fmt::print(stderr, "Gurobi Exception {}: {}\n", ex.getErrorCode(), ex.getMessage());
        throw ex;
    }
}

SolveState solveDominatingSet(PdsState& state, bool output, double timeLimit) {
    auto model = GRBModel(getEnv());
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    map <PowerGrid::VertexDescriptor, GRBVar> xi;
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

    return solveModel(state, xi, model);
}

SolveState solveBrimkovExpanded(PdsState& state, bool output, double timeLimit) {
    auto model = GRBModel(getEnv());
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    map <PowerGrid::VertexDescriptor, GRBVar> xi;
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

    return solveModel(state, xi, model);
}

SolveState solveBrimkov(PdsState& state, bool output, double timeLimit) {
    auto model = GRBModel(getEnv());
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    map <PowerGrid::VertexDescriptor, GRBVar> xi;
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

    return solveModel(state, xi, model);
}

SolveState solveJovanovic(PdsState& state, bool output, double timeLimit) {
    unused(state, output, timeLimit);
    struct NotImplemented : public std::exception {};
    throw NotImplemented{};
}

} // namespace pds
