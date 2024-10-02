#ifndef PDS_PATH_SOLVE
#define PDS_PATH_SOLVE

#include "pds.hpp"
#include "pdssolve.hpp"
#include "gurobi_common.hpp"

namespace pds {
namespace pathcallback {
enum class When {
    INTERMEDIATE_HS,
    FINAL
};
using PathCallback = std::function<void(When when, const PdsState& state, const std::vector<VertexList>& paths, size_t lower, size_t upper)>;
}
SolveResult solvePaths(PdsState &state,
                         int output,
                         double timeLimit,
                         int pathGenerator,
                         int pathInit,
                         int greedyUpper,
                         int earlyStop,
                         pathcallback::PathCallback pathCallback,
                         BoundCallback boundCallback,
                         int intermediatePaths,
                         int pathLimit);
} // namespace pds
#endif

