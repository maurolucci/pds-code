//
// Created by max on 19.08.22.
//

#ifndef PDS_PDSSOLVE_HPP
#define PDS_PDSSOLVE_HPP

#include "pds.hpp"
#include <concepts>

namespace pds {
template<std::invocable<const PdsState&, const std::string&> F>
bool exhaustiveSimpleReductions(PdsState& state, F callback = pds::unused) {
    bool anyChanged = false;
    bool changed;
    do {
        changed = false;
        if (state.disableLowDegree()) { callback(state, "low_degree"); changed = true; }
        while (state.collapseLeaves()) { callback(state, "leaves"); changed = true; }
        if (state.reduceObservedNonZi()) { callback(state, "non_zi"); changed = true; }
        while (state.collapseDegreeTwo()) { callback(state, "path"); changed = true; }
        if (state.collapseObservedEdges()) { callback(state, "observed_edges"); changed = true; }
        anyChanged |= changed;
    } while (changed);
    return anyChanged;
}

template<std::invocable<const PdsState&, const std::string&> F>
bool exhaustiveReductions(PdsState& state, F callback = unused) {
    bool anyChanged = false;
    bool changed;
    do {
        changed = exhaustiveSimpleReductions(state, callback);
        if (state.disableObservationNeighborhood()) { callback(state, "observation_neighborhood"); changed = true; }
        if (state.activateNecessaryNodes()) { callback(state, "necessary_nodes"); changed = true; }
    } while (changed);
    return anyChanged;
}

template<std::invocable<const PdsState&, const std::string&> F = void(const PdsState&, const std::string&)>
void solveGreedy(PdsState& state, F callback = unused) {
    using Vertex = PdsState::Vertex;
    while (!state.allObserved()) {
        exhaustiveReductions(state, callback);
        if (state.allObserved()) break;
        std::optional<Vertex> best;
        size_t bestObserved = 0;
        set<Vertex> active;
        for (auto v: state.graph().vertices()) {
            if (state.isActive(v)) { active.insert(v); }
        }
        for (auto v: state.graph().vertices()) {
            if (state.isBlank(v)) {
                size_t numObserved = observationNeighborhood(state.graph(), active).size();
                if (!best || bestObserved < numObserved) {
                    best = { v };
                    bestObserved = numObserved;
                }
            }
        }
        if (!best) break;
        state.setActive(*best);
    }
    return state;
}
}

#endif //PDS_PDSSOLVE_HPP
