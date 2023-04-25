#ifndef PDS_FORT_SOLVE
#define PDS_FORT_SOLVE

#include "pds.hpp"

namespace pds {
namespace callback {
enum class When {
    INTERMEDIATE_HS,
    FINAL
};
using FortCallback = std::function<void(When when, const PdsState& state, const std::vector<VertexList>& forts, size_t lower, size_t upper)>;
}
SolveResult solveBozeman(
        PdsState &state,
        int output,
        double timeLimit,
        int fortGenerator,
        int fortInit,
        int greedyUpper,
        int earlyStop,
        callback::FortCallback callback = noop_v);
} // namespace pds
#endif

