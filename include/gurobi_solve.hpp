//
// Created by max on 25.07.22.
//

#ifndef PDS_GUROBI_SOLVE_HPP
#define PDS_GUROBI_SOLVE_HPP

#include <range/v3/all.hpp>

#include "pds.hpp"

namespace pds {

SolveState solve_pds(PdsState&, bool output = false, double timeLimit = 10 * 60);

SolveState solveDominatingSet(PdsState& state, bool output = false, double timeLimit = 10 * 60);

SolveState solveBrimkovExpanded(PdsState& state, bool output = false, double timeLimit = 10 * 60);

SolveState solveBrimkov(PdsState& state, bool output = false, double timeLimit = 10 * 60);

SolveState solveJovanovic(PdsState& state, bool output = false, double timeLimit = 10 * 60);
} //namespace pds
#endif //PDS_GUROBI_SOLVE_HPP
