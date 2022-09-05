//
// Created by max on 10.08.22.
//
#include "gurobi_solve.hpp"

#include <gurobi_c++.h>

namespace pds {
bool solve_pds(PdsState& state, bool output, double timeLimit) {
    auto active = state.active();
    auto observed = state.observed();
    bool result = solve_pds(state.graph(), active, output, timeLimit);
    for (auto v: state.graph().vertices()) {
        if (active[v] == PmuState::Active) {
            state.setActive(v);
        }
    }
    return result;
}

bool solve_pds(const PowerGrid &graph, map <PowerGrid::vertex_descriptor, PmuState> &active, bool output, double timeLimit) {
    namespace r3 = ranges;
    auto env = GRBEnv();
    auto model = GRBModel(env);
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    //GRBVar *xi_p = model.addVars(num_vertices(graph), GRB_BINARY);
    //GRBVar *pij_p = model.addVars(2 * num_edges(graph), GRB_BINARY);
    //GRBVar *si_p = model.addVars(num_vertices(graph));
    map <PowerGrid::vertex_descriptor, GRBVar> xi;
    map <PowerGrid::vertex_descriptor, GRBVar> si;
    map <std::pair<PowerGrid::vertex_descriptor, PowerGrid::vertex_descriptor>, GRBVar> pij;
    const double M = 2 * graph.numVertices();
    {
        GRBLinExpr objective{};
        for (auto v: graph.vertices()) {
            xi.insert({v, model.addVar(0.0, M, 1.0, GRB_BINARY)});
            // si[v] >= 1 && si[v] <= num_vertices(graph)
            si.insert({v, model.addVar(1.0, M, 0.0, GRB_CONTINUOUS)});
            objective += xi[v];
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
            switch (active.at(w)) {
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
        if (active[v] != PmuState::Active) {
            model.addConstr(si[v] <= M * observingNeighbors);
        }
        model.addConstr(inObserved <= 1);
        if (graph[v].zero_injection) {
            model.addConstr(outObserved <= 1);
        } else {
            model.addConstr(outObserved == 0);
        }

        model.addConstr(si[v] <= graph.numVertices());

        if (active[v] == PmuState::Active) {
            model.addConstr(xi[v] == 1);
        } else if (active[v] == PmuState::Inactive) {
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

    model.optimize();
    if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
        return false;
    } else {
        std::cout << "result: " << model.get(GRB_DoubleAttr_ObjBound) << std::endl;
        for (auto v: graph.vertices()) {
            if (xi.at(v).get(GRB_DoubleAttr_X) > 0.5) {
                active[v] = PmuState::Active;
            } else {
                active[v] = PmuState::Inactive;
            }
        }
        return true;
    }
}

bool solveDominatingSet(PdsState& state, bool output, double timeLimit) {
    auto env = GRBEnv();
    auto model = GRBModel(env);
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    map <PowerGrid::vertex_descriptor, GRBVar> xi;
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

    model.optimize();
    if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
        return false;
    } else {
        std::cout << "result: " << model.get(GRB_DoubleAttr_ObjBound) << std::endl;
        for (auto v: state.graph().vertices()) {
            if (xi.at(v).get(GRB_DoubleAttr_X) > 0.5) {
                state.setActive(v);
            } else {
                state.setInactive(v);
            }
        }
    }
    return true;
}

bool solveBrimkov(PdsState& state, bool output, double timeLimit) {
    auto env = GRBEnv();
    auto model = GRBModel(env);
    model.set(GRB_IntParam_LogToConsole, int{output});
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    map <PowerGrid::vertex_descriptor, GRBVar> xi;
    map <PowerGrid::vertex_descriptor, GRBVar> si;
    map <PowerGrid::edge_descriptor, GRBVar> ye;
    GRBLinExpr objective;
    size_t T = state.graph().numVertices();
    for (auto v: state.graph().vertices()) {
        xi.try_emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY));
        si.try_emplace(v, model.addVar(0.0, static_cast<double>(T), 0.0, GRB_INTEGER));
        objective += xi.at(v);
        for (auto e: state.graph().outEdges(v)) {
            assert(!ye.contains(e));
            ye.try_emplace(e, model.addVar(0.0, 1.0, 0.0, GRB_BINARY));
        }
    }
    model.setObjective(objective, GRB_MINIMIZE);
    // x_v + \sum_{e \in N(v)} y_e = 1 (3)
    for (auto v: state.graph().vertices()) {
        GRBLinExpr maxPropagation = xi.at(v);
        for (auto e: state.graph().inEdges(v)) {
            unused(e);
            maxPropagation += ye.at(e);
        }
        model.addConstr(maxPropagation == 1);
    }
    // s_u - s_v + (T + 1) y_{uv} <= T (4) e@(u,v) \in E'
    // s_w - s_v + (T + 1) y_e <= T + (T+1) x_u , e@(u,v) \in E', w \in N(u) - v
    for (auto u: state.graph().vertices()) {
        for (auto e: state.graph().outEdges(u)) {
            if (!state.isZeroInjection(u)) model.addConstr(ye.at(e) == 0);
            auto v = state.graph().target(e);
            model.addConstr(ye.at(e) <= xi.at(u));
            //model.addConstr(si.at(u) - si.at(v) + (T + 1) * ye.at(e) <= T);
            for (auto w: state.graph().neighbors(u)) {
                if (u != v && u != w) {
                    unused(u, w, v, e);
                    //model.addConstr(si.at(w) - si.at(v) + (T + 1) * ye.at(e) <= T + (T+1) * xi.at(u));
                }
            }
        }
    }

    model.optimize();
    if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
        return false;
    } else {
        std::cout << "result: " << model.get(GRB_DoubleAttr_ObjBound) << std::endl;
        for (auto v: state.graph().vertices()) {
            if (xi.at(v).get(GRB_DoubleAttr_X) > 0.5) {
                state.setActive(v);
            } else {
                state.setInactive(v);
            }
        }
    }
    return true;
}

bool solveJovanovic(PdsState& state, bool output, double timeLimit) {
    unused(state, output, timeLimit);
    struct NotImplemented : public std::exception {};
    throw NotImplemented{};
}

} // namespace pds
