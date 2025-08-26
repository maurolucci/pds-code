#ifndef PDS_EFPS_SOLVE
#define PDS_EFPS_SOLVE

#include "gurobi_common.hpp"
#include "pds.hpp"
#include "pdssolve.hpp"

namespace pds {
namespace efpscallback {
enum class When { INTERMEDIATE_HS, FINAL };
using EFPSCallback =
    std::function<void(When when, const PdsState &state,
                       const std::vector<std::pair<EdgeList, size_t>> &efpss,
                       size_t lower, size_t upper)>;
} // namespace efpscallback
SolveResult solveEFPS(PdsState &state, int output, double timeLimit,
                      int efpsGenerator, int efpsInit, int greedyUpper,
                      int earlyStop, efpscallback::EFPSCallback efpsCallback,
                      BoundCallback boundCallback, int intermediateEFPS,
                      int efpsLimit);
} // namespace pds
#endif
