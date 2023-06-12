//
// Created by max on 25.07.22.
//

#ifndef PDS_GUROBI_SOLVE_HPP
#define PDS_GUROBI_SOLVE_HPP

#include <range/v3/all.hpp>

#include "pds.hpp"
#include "pdssolve.hpp"
#include "gurobi_common.hpp"

// Forward declaration to avoid full gurobi header inclusion
namespace pds {

inline const double TIME_LIMIT = 10 * 60;

MIPModel modelJovanovic(PdsState& state);
MIPModel modelJovanovicExpanded(PdsState& state);
MIPModel modelDomination(PdsState& state);
MIPModel modelBrimkov(PdsState& state);
MIPModel modelBrimkovExpanded(PdsState& state);
MIPModel modelAzamiBrimkov(PdsState& state);

SolveResult solveMIP(const PdsState& state, MIPModel&, bool output = false, double timeLimit = TIME_LIMIT, pds::BoundCallback = noop_v, bool validate = false);

void applySolution(PdsState&, MIPModel& model);

template<class F = decltype(modelJovanovicExpanded)>
//requires std::is_invocable_v<PdsState&>
inline SolveResult solvePowerDominatingSet(PdsState& state, bool output, double timeLimit, pds::BoundCallback boundCallback = noop_v, F model = modelJovanovicExpanded, bool validate = false) {
    auto mip = model(state);
    auto result = solveMIP(state, mip, output, timeLimit, boundCallback, validate);
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
