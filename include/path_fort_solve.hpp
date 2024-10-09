#ifndef PDS_PATH_FORT_SOLVE
#define PDS_PATH_FORT_SOLVE

#include "pds.hpp"
#include "path_solve.hpp"
#include "fort_solve.hpp"
#include "pdssolve.hpp"
#include "gurobi_common.hpp"

namespace pds {
SolveResult solvePathsForts(
                         PdsState &state,
                         int output,
                         double timeLimit,
                         int pathGenerator,
                         int pathInit,
                         int greedyUpper,
                         int earlyStop,
                         callback::FortCallback fortCallback,
                         pathcallback::PathCallback pathCallback,
                         BoundCallback boundCallback,
                         int intermediatePaths,
                         int pathLimit);
} // namespace pds
#endif

