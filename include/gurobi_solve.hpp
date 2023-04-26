//
// Created by max on 25.07.22.
//

#ifndef PDS_GUROBI_SOLVE_HPP
#define PDS_GUROBI_SOLVE_HPP

#include <range/v3/all.hpp>

#include "pds.hpp"
#include "pdssolve.hpp"

// Forward declaration to avoid full gurobi header inclusion
struct GRBModel;
struct GRBVar;
struct GRBEnv;

namespace pds {

inline const double TIME_LIMIT = 10 * 60;

struct MIPModel {
    std::unique_ptr<GRBModel> model;
    map<PowerGrid::VertexDescriptor, GRBVar> xi;
    MIPModel();
    MIPModel(MIPModel&& other) = default;
    virtual ~MIPModel();
};

void preloadMIPSolver();
GRBEnv& getEnv();

void relaxMIPModel(MIPModel&);

MIPModel modelJovanovic(PdsState& state);
MIPModel modelJovanovicExpanded(PdsState& state);
MIPModel modelDomination(PdsState& state);
MIPModel modelBrimkov(PdsState& state);
MIPModel modelBrimkovExpanded(PdsState& state);
MIPModel modelAzamiBrimkov(PdsState& state);

SolveResult solveMIP(MIPModel&, bool output = false, double timeLimit = TIME_LIMIT, pds::BoundCallback = noop_v);

void applySolution(PdsState&, MIPModel& model);

template<class F = decltype(modelJovanovicExpanded)>
//requires std::is_invocable_v<PdsState&>
inline SolveResult solvePowerDominatingSet(PdsState& state, bool output, double timeLimit, pds::BoundCallback boundCallback = noop_v, F model = modelJovanovicExpanded) {
    auto mip = model(state);
    auto result = solveMIP(mip, output, timeLimit, boundCallback);
    if (result.state != SolveState::Infeasible) {
        applySolution(state, mip);
    }
    return result;
}

//inline SolveState solve_pds(PdsState& state, bool output = false, double timeLimit = TIME_LIMIT) {
//    return solve_pds(state, output, timeLimit, modelJovanovicExpanded);
//}

} //namespace pds
#endif //PDS_GUROBI_SOLVE_HPP
