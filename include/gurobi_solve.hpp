//
// Created by max on 25.07.22.
//

#ifndef PDS_GUROBI_SOLVE_HPP
#define PDS_GUROBI_SOLVE_HPP

#include <range/v3/all.hpp>

#include "pds.hpp"

namespace pds {

inline const double TIME_LIMIT = 10 * 60;

/// Opaque implementation does not leak internal usage of Gurobi
struct MIPModel {
    struct Implementation;
    std::unique_ptr<pds::MIPModel::Implementation> impl;
    MIPModel() = default;
    MIPModel(MIPModel&& other) = default;
    virtual ~MIPModel();
};

void preloadMIPSolver();

void relaxMIPModel(MIPModel&);

MIPModel modelJovanovic(PdsState& state);
MIPModel modelJovanovicExpanded(PdsState& state);
MIPModel modelDomination(PdsState& state);
MIPModel modelBrimkov(PdsState& state);
MIPModel modelBrimkovExpanded(PdsState& state);
MIPModel modelAzamiBrimkov(PdsState& state);

SolveState solveMIP(MIPModel&, bool output = false, double timeLimit = TIME_LIMIT);

void applySolution(PdsState&, MIPModel& model);

template<class F = decltype(modelJovanovicExpanded)>
//requires std::is_invocable_v<PdsState&>
inline SolveState solvePowerDominatingSet(PdsState& state, bool output, double timeLimit, F model = modelJovanovicExpanded) {
    auto mip = model(state);
    auto result = solveMIP(mip, output, timeLimit);
    if (result != SolveState::Infeasible) {
        applySolution(state, mip);
    }
    return result;
}

//inline SolveState solve_pds(PdsState& state, bool output = false, double timeLimit = TIME_LIMIT) {
//    return solve_pds(state, output, timeLimit, modelJovanovicExpanded);
//}

} //namespace pds
#endif //PDS_GUROBI_SOLVE_HPP
