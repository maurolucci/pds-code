//
// Created by max on 25.07.22.
//

#ifndef PDS_GUROBI_SOLVE_HPP
#define PDS_GUROBI_SOLVE_HPP

#include <gurobi_c++.h>
#include <range/v3/all.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "utility.hpp"
#include "pds.hpp"

namespace pds {

template<class Graph, class ActivePropertyMap, class ObservedPropertyMap>
bool block_deg2(Graph& graph, ActivePropertyMap& active, ObservedPropertyMap& observed) {
    bool has_deg3 = ranges::any_of(
            range_pair(vertices(graph)),
            [&graph](const auto& v) { return degree(v, graph) > 3; });
    for (const auto& v: range_pair(vertices(graph))) {
        if (has_deg3 && boost::degree(v, graph) == 2) {
            active[v] = PmuState::Inactive;
        }
    }
}

template<class Graph, class ActivePropertyMap, class ObservedPropertyMap>
bool solve_pds(const Graph &graph, ActivePropertyMap &active, ObservedPropertyMap &observed) {
    namespace r3 = ranges;
    auto env = GRBEnv();
    auto model = GRBModel(env);
    //GRBVar *xi_p = model.addVars(num_vertices(graph), GRB_BINARY);
    //GRBVar *pij_p = model.addVars(2 * num_edges(graph), GRB_BINARY);
    //GRBVar *si_p = model.addVars(num_vertices(graph));
    std::map<typename Graph::vertex_descriptor, GRBVar> xi;
    std::map<typename Graph::vertex_descriptor, GRBVar> si;
    std::map<std::pair<typename Graph::vertex_descriptor, typename Graph::vertex_descriptor>, GRBVar> pij;
    const double M = 2 * num_vertices(graph);
    {
        GRBLinExpr objective{};
        for (auto [i, v]: range_pair(vertices(graph)) | r3::views::enumerate) {
            xi[v] = model.addVar(0.0, M, 1.0, GRB_BINARY);
            // si[v] >= 1 && si[v] <= num_vertices(graph)
            si[v] = model.addVar(1.0, M, 0.0, GRB_CONTINUOUS);
            objective += xi[v];
        }
        for (auto [i, e]: range_pair(edges(graph)) | r3::views::enumerate) {
            auto v = source(e, graph);
            auto w = target(e, graph);
            pij[{v, w}] = model.addVar(0.0, M, 0.0, GRB_BINARY);//pij_p[2 * i];
            pij[{w, v}] = model.addVar(0.0, M, 0.0, GRB_BINARY);//pij_p[2 * i + 1];
        }
        //model.setObjective(objective, GRB_MINIMIZE);
    }

    for (auto v: range_pair(vertices(graph))) {
        model.addConstr(si[v] >= 1);
        for (auto w: r3::concat_view(range_pair(adjacent_vertices(v, graph)), r3::single_view{v})) {
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
        for (auto e: range_pair(out_edges(v, graph))) {
            auto w = target(e, graph);
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

        model.addConstr(si[v] <= num_vertices(graph));

        if (active[v] == PmuState::Active) {
            model.addConstr(xi[v] == 1);
        } else if (active[v] == PmuState::Inactive) {
            model.addConstr(xi[v] == 0);
        }
    }

    for (auto e: range_pair(edges(graph))) {
        auto v = source(e, graph);
        auto w = target(e, graph);
        model.addConstr(pij[{v, w}] + pij[{w, v}] <= 1);
        for (auto t: r3::concat_view(range_pair(adjacent_vertices(w, graph)), r3::single_view{w})) {
            if (t != v)
                model.addConstr(si[v] >= si[t] + 1 - M * (1 - pij[{w, v}]));
        }
        for (auto t: r3::concat_view(range_pair(adjacent_vertices(v, graph)), r3::single_view{v})) {
            if (t != w)
                model.addConstr(si[w] >= si[t] + 1 - M * (1 - pij[{v, w}]));
        }
    }

    model.optimize();
    std::cout << "result: " << model.get(GRB_DoubleAttr_ObjBound) << std::endl;
    for (auto v: range_pair(vertices(graph))) {
        if (xi[v].get(GRB_DoubleAttr_X) > 0.5) {
            active[v] = PmuState::Active;
        }
    }
    return false;
}

} //namespace pds
#endif //PDS_GUROBI_SOLVE_HPP
