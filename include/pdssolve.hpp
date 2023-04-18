//
// Created by max on 19.08.22.
//

#ifndef PDS_PDSSOLVE_HPP
#define PDS_PDSSOLVE_HPP

#include "pds.hpp"
#include <queue>
#include <concepts>

namespace pds {
template<std::invocable<const PdsState&, const std::string&> F = void(const PdsState&, const std::string&)>
bool exhaustiveSimpleReductions(PdsState& state, F callback = pds::unused) {
    bool anyChanged = false;
    bool changed;
    do {
        changed = false;
        if (state.disableLowDegree()) { callback(state, "low_degree"); changed = true; }
        while (state.collapseLeaves()) { callback(state, "leaves"); changed = true; }
        while (state.collapseDegreeTwo()) { callback(state, "path"); changed = true; }
        if (state.reduceObservedNonZi()) { callback(state, "non_zi"); changed = true; }
        if (state.collapseObservedEdges()) { callback(state, "observed_edges"); changed = true; }
        anyChanged |= changed;
    } while (changed);
    return anyChanged;
}

template<std::invocable<const PdsState&, const std::string&> F = void(const PdsState&, const std::string&)>
bool dominationReductions(PdsState &state, bool firstRun = true, F callback = unused) {
    bool changed = false;
    if (state.disableObservationNeighborhood()) { callback(state, "observation_neighborhood"); changed = true; }
    if ((firstRun || changed) && state.activateNecessaryNodes()) { callback(state, "necessary_nodes"); changed = true; }
    return changed;
}

template<std::invocable<const PdsState&, const std::string&> F = void(const PdsState&, const std::string&)>
bool exhaustiveReductions(PdsState &state, bool firstRun = true, F callback = unused) {
    bool anyChanged = false;
    bool changed;
    do {
        changed = exhaustiveSimpleReductions(state, callback);
        if (firstRun || changed) changed |= dominationReductions(state, firstRun, callback);
        firstRun = false;
        anyChanged |= changed;
    } while (changed);
    return anyChanged;
}

template<std::invocable<const PdsState&, const std::string&> F = void(const PdsState&, const std::string&)>
bool noNecessaryReductions(PdsState &state, bool firstRun = true, F callback = unused) {
    bool anyChanged = false;
    bool changed;
    do {
        changed = exhaustiveSimpleReductions(state, callback);
        if (firstRun || changed) if(state.disableObservationNeighborhood()) {
            callback(state, "observation_neighborhood");
            changed = true;
        }
        firstRun = false;
        anyChanged |= changed;
    } while (changed);
    return anyChanged;
}

namespace greedy_strategies {
std::optional<PdsState::Vertex> largestObservationNeighborhood(const PdsState &state);

std::optional<PdsState::Vertex> largestDegree(const PdsState &state);

std::optional<PdsState::Vertex> medianDegree(const PdsState &state);
}

template<
        std::invocable<const PdsState&> Strategy = std::optional<PdsState::Vertex>(const PdsState&),
        std::invocable<const PdsState&, const std::string&> F = void(const PdsState&, const std::string&)
>
SolveResult solveGreedy(PdsState& state, bool applyReductions = true, Strategy strategy = greedy_strategies::largestObservationNeighborhood, F callback = unused) {
    while (!state.allObserved()) {
        if (applyReductions) {
            exhaustiveReductions(state, true, callback);
        }
        if (state.allObserved()) break;
        auto best = strategy(state);
        if (!best) break;
        state.setActive(*best);
    }
    return {size_t{0}, state.numActive(), SolveState::Heuristic};
}

SolveResult fastGreedy(PdsState& state);

SolveResult topDownGreedy(PdsState &state);

using Bounds = std::pair<size_t, size_t>;
Bounds sensorBounds(const PdsState& state);

template< std::invocable<const PdsState&> Strategy = std::optional<PdsState::Vertex>(const PdsState&) >
SolveResult solveBranching(PdsState &state,
                    bool useReductions,
                    Strategy strategy = greedy_strategies::largestDegree) {
    if (useReductions) {
        exhaustiveReductions(state, true);
    }
    auto heuristic = state;
    fastGreedy(heuristic);
    auto upper = sensorBounds(heuristic).second;
    fmt::print("heuristic result: {}\n", upper);
    size_t lower = 0;
    auto compare = [](const auto& first, const auto& second) { return first.first.first > second.first.first; };
    using Element = std::pair<Bounds, PowerGrid>;
    std::priority_queue<Element, std::vector<Element>, decltype(compare)> queue(compare);
    queue.push({sensorBounds(state), state.graph()});
    size_t explored = 0;
    using namespace std::chrono_literals;
    auto now = []() { return std::chrono::high_resolution_clock::now(); };
    auto sec = [](auto time) { return std::chrono::duration_cast<std::chrono::seconds>(time).count(); };
    auto printPeriod = 1s;
    auto previousPrint = now();
    auto start = previousPrint;
    while (!queue.empty()) {
        ++explored;
        PdsState top(std::move(queue.top().second));
        auto bounds = queue.top().first;
        queue.pop();
        if (bounds.first > upper) continue;
        lower = bounds.first;
        auto t = now();
        if (t - previousPrint > printPeriod) {
            fmt::print("explored {} nodes\t{}\t{}\t{}\t{}\t{}\t{}s\n", explored, lower, upper, bounds.first, bounds.second, top.allObserved(), sec(t - start));
            previousPrint = t;
        }
        upper = std::min(upper, bounds.second);
        if (bounds.first == bounds.second && top.allObserved()) {
            heuristic = top;
            fmt::print("incumbent solution: {}\t{}\n", bounds.first, top.allObserved());
        }
        auto best = strategy(top);
        if (!best) continue;
        auto activated = top;
        activated.setActive(*best);
        top.setInactive(*best);
        if (useReductions) {
            exhaustiveReductions(activated);
            exhaustiveReductions(top);
        }
        auto activatedBounds = sensorBounds(activated);
        auto disabledBounds = sensorBounds(top);
        if (activatedBounds.first < upper) {// && isFeasible(activated)) {// && activatedBounds.first <= activatedBounds.second
            upper = std::min(upper, activatedBounds.second);
            queue.template emplace(activatedBounds, std::move(activated.moveGraph()));
        }
        if (disabledBounds.first < upper) {// && isFeasible(top)) { // && disabledBounds.first <= disabledBounds.second
            upper = std::min(upper, disabledBounds.second);
            queue.emplace(disabledBounds, std::move(top.moveGraph()));
        }
    }
    fmt::print("finished after exploring {} nodes\t{}\t{}\n", explored, lower, upper);
    state = std::move(heuristic);
    fmt::print("solved by branching. result: {}\n", upper);
    return {lower, upper, SolveState::Optimal};
}

} //namespace pds

#endif //PDS_PDSSOLVE_HPP
