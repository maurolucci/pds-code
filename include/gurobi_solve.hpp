//
// Created by max on 25.07.22.
//

#ifndef PDS_GUROBI_SOLVE_HPP
#define PDS_GUROBI_SOLVE_HPP

#include <gurobi_c++.h>
#include <range/v3/all.hpp>

#include "pds.hpp"

namespace pds {

inline bool solve_pds(const PowerGrid &graph, map<PowerGrid::vertex_descriptor, PmuState> &active) {
    namespace r3 = ranges;
    auto env = GRBEnv();
    auto model = GRBModel(env);
    model.set(GRB_IntParam_LogToConsole, 0);
    model.set(GRB_StringParam_LogFile, "gurobi.log");
    //GRBVar *xi_p = model.addVars(num_vertices(graph), GRB_BINARY);
    //GRBVar *pij_p = model.addVars(2 * num_edges(graph), GRB_BINARY);
    //GRBVar *si_p = model.addVars(num_vertices(graph));
    map<PowerGrid::vertex_descriptor, GRBVar> xi;
    map<PowerGrid::vertex_descriptor, GRBVar> si;
    map<std::pair<PowerGrid::vertex_descriptor, PowerGrid::vertex_descriptor>, GRBVar> pij;
    const double M = 2 * graph.numVertices();
    {
        GRBLinExpr objective{};
        for (auto v: graph.vertices()) {
            xi[v] = model.addVar(0.0, M, 1.0, GRB_BINARY);
            // si[v] >= 1 && si[v] <= num_vertices(graph)
            si[v] = model.addVar(1.0, M, 0.0, GRB_CONTINUOUS);
            objective += xi[v];
        }
        for (auto e: graph.edges()) {
            auto v = graph.source(e);
            auto w = graph.target(e);
            pij[{v, w}] = model.addVar(0.0, M, 0.0, GRB_BINARY);//pij_p[2 * i];
            pij[{w, v}] = model.addVar(0.0, M, 0.0, GRB_BINARY);//pij_p[2 * i + 1];
        }
        //model.setObjective(objective, GRB_MINIMIZE);
    }

    for (auto v: graph.vertices()) {
        model.addConstr(si[v] >= 1);
        for (auto w: r3::concat_view(graph.neighbors(v), r3::single_view{v})) {
            switch (active[v]) {
            case PmuState::Blank:
                model.addConstr(si[v] <= xi[w] + M * (1 - xi[w]));
                break;
            case PmuState::Inactive:
                model.addConstr(si[v] <= M);
                break;
            case PmuState::Active:
                model.addConstr(si[v] <= 1);
            }
        }
        GRBLinExpr observingNeighbors = xi[v];
        // v observes at most one neighbor w1
        GRBLinExpr outObserved;
        // v is observed by at most one neighbor w1
        GRBLinExpr inObserved;
        for (auto w: graph.neighbors(v)) {
            observingNeighbors += xi[w];
            if (graph[w].zero_injection) {
                observingNeighbors += pij[{w, v}];
            }

            outObserved += pij[{v, w}];
            inObserved += pij[{w, v}];
        }
        if (active[v] != PmuState::Active) {
            model.addConstr(si[v] <= M * observingNeighbors);
        }
        model.addConstr(inObserved <= 1);
        model.addConstr(outObserved <= 1);

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
        model.addConstr(pij[{v, w}] + pij[{w, v}] <= 1);
        for (auto t: r3::concat_view(graph.neighbors(w), r3::single_view{w})) {
            if (t != v)
                model.addConstr(si[v] >= si[t] + 1 - M * (1 - pij[{w, v}]));
        }
        for (auto t: r3::concat_view(graph.neighbors(v), r3::single_view{v})) {
            if (t != w)
                model.addConstr(si[w] >= si[t] + 1 - M * (1 - pij[{v, w}]));
        }
    }

    model.optimize();
    std::cout << "result: " << model.get(GRB_DoubleAttr_ObjBound) << std::endl;
    for (auto v: graph.vertices()) {
        if (xi[v].get(GRB_DoubleAttr_X) > 0.5) {
            active[v] = PmuState::Active;
        } else {
            active[v] = PmuState::Inactive;
        }
    }
    return false;
}

} //namespace pds
#endif //PDS_GUROBI_SOLVE_HPP
