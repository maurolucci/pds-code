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

SolveResult fastGreedy(PdsState &state) {
    //auto comp = [&state] (auto v, auto w) { return state.graph().degree(v) < state.graph().degree(w); };
    auto vertices = state.graph().vertices()
                    | ranges::views::filter([&state](auto v) { return state.isBlank(v); })
                    | ranges::views::transform([&state](auto v) { return std::pair(state.graph().degree(v), v); })
                    | ranges::to<std::vector>;
    ranges::make_heap(vertices);
    while (!state.allObserved() && !vertices.empty()) {
        while (!state.graph().hasVertex(vertices.front().second) || !state.isBlank(vertices.front().second)) {
            ranges::pop_heap(vertices);
            vertices.pop_back();
            if (vertices.empty()) break;
        }
        ranges::pop_heap(vertices);
        auto v = vertices.back().second;
        vertices.pop_back();
        state.setActive(v);
        exhaustiveReductions(state);
    }
    if (!state.allObserved()) return { state.graph().numVertices(), size_t{0}, SolveState::Infeasible };
    return {size_t{0}, state.numObserved(), SolveState::Heuristic};
}

SolveResult topDownGreedy(PdsState &state) {
    auto vertices = state.graph().vertices()
                    | ranges::views::filter([&state](auto v) { return state.isBlank(v); })
                    | ranges::views::transform([&state](auto v) { return std::pair(-ssize_t(state.graph().degree(v)), v); })
                    | ranges::to<std::vector>;
    for (auto v: vertices) {
        state.setActive(v.second);
    }
    if (!state.allObserved()) return { state.graph().numVertices(), size_t{0}, SolveState::Infeasible };
    ranges::make_heap(vertices);
    while (!vertices.empty()) {
        ranges::pop_heap(vertices);
        auto v = vertices.back().second;
        vertices.pop_back();
        state.setInactive(v);
        if (!state.allObserved()) {
            state.setActive(v);
        }
    }
    return {size_t{0}, state.numObserved(), SolveState::Heuristic};
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