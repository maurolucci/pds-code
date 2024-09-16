#ifndef PDS_FORT_SOLVE
#define PDS_FORT_SOLVE

#include "pds.hpp"
#include "pdssolve.hpp"
#include "gurobi_common.hpp"

namespace pds {
namespace callback {
enum class When {
    INTERMEDIATE_HS,
    FINAL
};
using CycleCallback = std::function<void(When when, const PdsState& state, const std::vector<VertexList>& cycles, size_t lower, size_t upper)>;
}
SolveResult solveCycles(PdsState &state,
                         int output,
                         double timeLimit,
                         int cycleGenerator,
                         int cycleInit,
                         int greedyUpper,
                         int earlyStop,
                         callback::CycleCallback cycleCallback,
                         BoundCallback boundCallback,
                         int intermediateCycles);
SolveResult solveLazyForts(PdsState& state, int output, double timeLimit, callback::CycleCallback CycleCB, BoundCallback boundsCB);
} // namespace pds
#endif

