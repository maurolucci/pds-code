//
// Created by max on 06.02.23.
//

#include "pdssolve.hpp"
namespace pds {


Bounds sensorBounds(const PdsState &state) {
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

SolveResult fastGreedy(PdsState &state, int useReductions) {
    if (state.allObserved()) { return {size_t{0}, state.numObserved(), SolveState::Heuristic}; }
    //auto comp = [&state] (auto v, auto w) { return state.graph().degree(v) < state.graph().degree(w); };
    if (state.allObserved()) {
        return {state.numActive(), state.numActive(), SolveState::Optimal};
    }
    size_t lower = state.numActive();
    auto blank = state.graph().vertices()
                    | ranges::views::filter([&state](auto v) { return state.isBlank(v); })
                    | ranges::views::transform([&state](auto v) { return std::pair(state.unobservedDegree(v), v); })
                    | ranges::to<std::vector>;
    std::span vertices(blank);
    ranges::make_heap(vertices);
    while (!state.allObserved() && !vertices.empty()) {
        while (!state.graph().hasVertex(vertices.front().second) || !state.isBlank(vertices.front().second)) {
            ranges::pop_heap(vertices);
            vertices = vertices.subspan(0, vertices.size() - 1);
            if (vertices.empty()) break;
        }
        if (!vertices.empty()) {
            ranges::pop_heap(vertices);
            auto v = vertices.back().second;
            vertices = vertices.subspan(0, vertices.size() - 1);
            if (!state.isObserved(v) || state.unobservedDegree(v) != 0) {
                state.setActive(v);
                if (useReductions == 1) { exhaustiveReductions(state); }
                if (useReductions == 2) { dominationReductions(state); }
            }
        }
    }
    //size_t initial = state.numActive();
    for (auto [deg, v]: blank) {
        if (state.isActive(v)) {
            state.setBlank(v);
            if (!state.allObserved()) { state.setActive(v); }
        }
    }
    //fmt::print("greedy: {} → {}\n", initial, state.numActive());
    if (!state.allObserved()) return { state.graph().numVertices(), size_t{0}, SolveState::Infeasible };
    return {lower, state.numActive(), SolveState::Heuristic};
}

SolveResult topDownGreedy(PdsState &state, bool activateAll, std::span<PdsState::Vertex> vertices) {
    std::vector<PdsState::Vertex> blank;
    if (state.allObserved()) {
        return {state.numActive(), state.numActive(), SolveState::Optimal};
    }
    size_t lower = state.numActive() + 1;
    if (vertices.empty() && state.numBlank() > 0) {
        blank = state.graph().vertices()
                | ranges::views::filter([&state] (auto v) { return state.isBlank(v); })
                | ranges::to<std::vector<PdsState::Vertex>>;
        vertices = blank;
    }
    ranges::sort(vertices, [&state] (auto left, auto right) { return state.graph().degree(left) < state.graph().degree(right); });
    //auto vertices = state.graph().vertices()
    //                | ranges::views::filter([&state](auto v) { return state.isBlank(v); })
    //                | ranges::views::transform([&state](auto v) { return std::pair(-ssize_t(state.graph().degree(v)), v); })
    //                | ranges::to<std::vector>;
    if (activateAll) {
        for (auto v: vertices) {
            state.setActive(v);
        }
    }
    if (!state.allObserved()) return { state.graph().numVertices(), size_t{0}, SolveState::Infeasible };
    for (auto v: vertices) {
        if (state.isActive(v)) {
            state.setInactive(v);
            if (!state.allObserved()) {
                state.setActive(v);
            }
        }
    }
    return {lower, state.numActive(), SolveState::Heuristic};
}

namespace greedy_strategies {
std::optional<PdsState::Vertex> medianDegree(const PdsState &state) {
    using Vertex = PdsState::Vertex;
    std::vector<Vertex> vertices;
    for (auto v: state.graph().vertices()) {
        if (state.isBlank(v)) {
            vertices.push_back(v);
        }
    }
    if (vertices.size() == 0) return {};
    auto deg = [&state](auto v) { return state.graph().degree(v); };
    auto best = vertices.begin() + (vertices.size() + 1) / 2;
    ranges::nth_element(vertices, best, [deg](auto v, auto w) { return deg(v) < deg(w); });
    if (best != vertices.end()) {
        return {*best};
    } else {
        return {};
    }
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
            size_t numObserved = state.numObserved(); //TODO
            if (!best || bestObserved < numObserved) {
                best = {v};
                bestObserved = numObserved;
            }
        }
    }
    return best;
}

} // namespace greedy_strategies
} // namespace pds