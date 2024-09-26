#ifndef PDS_CYCLE_FORT_SOLVE
#define PDS_CYCLE_FORT_SOLVE

#include "pds.hpp"
#include "cycle_solve.hpp"
#include "fort_solve.hpp"
#include "pdssolve.hpp"
#include "gurobi_common.hpp"

namespace pds {
SolveResult solveCyclesForts(
                         PdsState &state,
                         int output,
                         double timeLimit,
                         int cycleGenerator,
                         int cycleInit,
                         int greedyUpper,
                         int earlyStop,
                         callback::FortCallback fortCallback,
                         cyclecallback::CycleCallback cycleCallback,
                         BoundCallback boundCallback,
                         int intermediateCycles);
} // namespace pds
#endif

