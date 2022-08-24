//
// Created by max on 19.08.22.
//

#ifndef PDS_PDSSOLVE_HPP
#define PDS_PDSSOLVE_HPP

#include "pds.hpp"
#include <queue>
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

template<std::invocable<const PdsState&, const std::string&> F = void(const PdsState&, const std::string&)>
bool exhaustiveReductions(PdsState &state, bool firstRun = true, F callback = unused) {
    bool anyChanged = false;
    bool changed;
    do {
        changed = exhaustiveSimpleReductions(state, callback);
        if ((firstRun || changed) && state.disableObservationNeighborhood()) { callback(state, "observation_neighborhood"); changed = true; }
        if ((firstRun || changed) && state.activateNecessaryNodes()) { callback(state, "necessary_nodes"); changed = true; }
        firstRun = false;
    } while (changed);
    return anyChanged;
}

namespace greedy_strategies {
std::optional<PdsState::Vertex> largestObservationNeighborhood(const PdsState &state) {
    using Vertex = PdsState::Vertex;
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
                best = {v};
                bestObserved = numObserved;
            }
        }
    }
    return best;
}

std::optional<PdsState::Vertex> largestDegree(const PdsState &state) {
    using Vertex = PdsState::Vertex;
    std::optional<Vertex> best;
    size_t bestObserved = 0;
    for (auto v: state.graph().vertices()) {
        if (state.isBlank(v)) {
            if (!best || bestObserved < state.graph().degree(v)) {
                best = {v};
                bestObserved = state.graph().degree(v);
            }
        }
    }
    return best;
}

std::optional<PdsState::Vertex> medianDegree(const PdsState &state) {
    using Vertex = PdsState::Vertex;
    std::vector<Vertex> vertices;
    for (auto v: state.graph().vertices()) {
        if (state.isBlank(v)) {
            vertices.push_back(v);
        }
    }
    if (vertices.size() == 0) return {};
    auto deg = [&state](auto v) { return state.graph().degree(v);};
    auto best = vertices.begin() + (vertices.size() + 1) / 2;
    ranges::nth_element(vertices, best, [deg](auto v, auto w) { return deg(v) < deg(w);});
    if (best != vertices.end()) {
        return {*best};
    } else {
        return {};
    }
}
}

template<
        std::invocable<const PdsState&> Strategy = std::optional<PdsState::Vertex>(const PdsState&),
        std::invocable<const PdsState&, const std::string&> F = void(const PdsState&, const std::string&)
>
void solveGreedy(PdsState& state, bool applyReductions = true, Strategy strategy = greedy_strategies::largestObservationNeighborhood, F callback = unused) {
    while (!state.allObserved()) {
        if (applyReductions) {
            exhaustiveReductions(state, true, callback);
        }
        if (state.allObserved()) break;
        auto best = strategy(state);
        if (!best) break;
        state.setActive(*best);
    }
}

using Bounds = std::pair<size_t, size_t>;
Bounds sensorBounds(const PdsState& state) {
    size_t lower = 0, upper = 0, unobserved = 0;
    for (auto v: state.graph().vertices()) {
        if (state.isActive(v)) {
            lower += 1;
            upper += 1;
        }
        if (!state.isObserved(v)) {
            unobserved += 1;
            if (state.isBlank(v)) {
                upper += 1;
            }
        }
    }
    return {lower + (unobserved > 0), upper + (unobserved > 0)};
}

bool isFeasible(const PdsState& state) {
    using Vertex = PdsState::Vertex;
    set<Vertex> active;
    for (auto v: state.graph().vertices()) {
        if (!state.isInactive(v)) { active.insert(v); }
    }
    return observationNeighborhood(state.graph(), active).size() == state.graph().numVertices();
}

template< std::invocable<const PdsState&> Strategy = std::optional<PdsState::Vertex>(const PdsState&) >
void solveBranching(PdsState &state,
                    bool useReductions,
                    Strategy strategy = greedy_strategies::largestDegree) {
    if (useReductions) {
        exhaustiveReductions(state, true);
    }
    auto heuristic = state;
    solveGreedy(heuristic, true, strategy);
    auto upper = sensorBounds(heuristic).second;
    fmt::print("heuristic result: {}\n", upper);
    size_t lower = 0;
    auto compare = [](const auto& first, const auto& second) { return 2 * first.first.first + first.first.second > 2 * second.first.first + second.first.second; };
    using Element = std::pair<Bounds, PdsState>;
    std::priority_queue<Element, std::vector<Element>, decltype(compare)> queue(compare);
    queue.push({sensorBounds(state), state});
    size_t explored = 0;
    while (!queue.empty()) {
        ++explored;
        auto [bounds, top] = queue.top();
        queue.pop();
        if (bounds.first > upper) continue;
        lower = bounds.first;
        fmt::print("explored {} nodes\t{}\t{}\t{}\t{}\t{}\n", explored, lower, upper, bounds.first, bounds.second, top.allObserved());
        upper = std::min(upper, bounds.second);
        if (bounds.first == bounds.second) {
            heuristic = top;
            fmt::print("incumbent solution: {}\t{}\n", bounds.first, top.allObserved());
        }
        auto best = strategy(top);
        if (!best) continue;
        auto activated = top;
        activated.setActive(*best);
        top.setInactive(*best);
        if (useReductions) {
            exhaustiveReductions(activated, true);
            exhaustiveReductions(top, true);
        }
        auto activatedBounds = sensorBounds(activated);
        auto disabledBounds = sensorBounds(top);
        if (activatedBounds.first < upper && isFeasible(activated)) {// && activatedBounds.first <= activatedBounds.second
            upper = std::min(upper, activatedBounds.second);
            queue.push({activatedBounds, activated});
        }
        if (disabledBounds.first < upper && isFeasible(top)) { // && disabledBounds.first <= disabledBounds.second
            upper = std::min(upper, disabledBounds.second);
            queue.push({disabledBounds, top});
        }
    }
    state = heuristic;
    fmt::print("solved by branching. result: {}\n", upper);
}

} //namespace pds

#endif //PDS_PDSSOLVE_HPP
