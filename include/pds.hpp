//
// Created by max on 19.07.22.
//

#ifndef PDS_PDS_HPP
#define PDS_PDS_HPP

#include <range/v3/all.hpp>

#include <iostream>
#include <fstream>
#include <optional>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include "map.hpp"
#include "utility.hpp"
#include <setgraph/graph.hpp>

namespace pds {

enum class PmuState : signed char {
    Blank = 0,
    Active = 1,
    Inactive = -1
};

struct Bus {
    std::string name;
    long id;
    bool zero_injection;
};

using PowerGrid = setgraph::SetGraph<Bus, setgraph::Empty, setgraph::EdgeDirection::Undirected>;

void exportGraphml(const PowerGrid& grid, std::ostream& out);

set<PowerGrid::vertex_descriptor> observationNeighborhood(const PowerGrid &graph, const set<PowerGrid::vertex_descriptor> &starts);

bool propagate(const PowerGrid&, set<PowerGrid::vertex_descriptor>&, size_t = 1);

set<PowerGrid::vertex_descriptor> observationNeighborhood(const PowerGrid& graph, const set<PowerGrid::vertex_descriptor>& starts);

template<class T>
size_t intersectionSize(const set<T>& first, const set<T>& second) {
    size_t count = 0;
    for (auto x: first) {
        if (second.contains(x)) {
            ++count;
        }
    }
    return count;
}

template<class T>
bool isSuperset(const set <T> &container, const set <T> &other) {
    if (other.size() > container.size()) return false;
    for (auto x: other) {
        if (!container.contains(x)) return false;
    }
    return true;
}

class PdsState {
public:
    using Vertex = PowerGrid::vertex_descriptor;
private:
    pds::map<Vertex, ssize_t> m_unobserved_degree;
    pds::set<Vertex> m_observed;
    pds::map<Vertex, PmuState> m_active;
    PowerGrid m_graph;

    std::vector<Vertex> m_steps_observed;
    std::vector<std::pair<Vertex, PmuState>> m_steps_pmu;
    std::vector<std::pair<size_t, size_t>> m_checkpoints;

    void propagate(Vertex vertex);

    void observe(Vertex vertex);

    inline bool disableLowDegreeRecursive(
            PowerGrid::vertex_descriptor start,
            set<PowerGrid::vertex_descriptor>& seen
    );

public:
    PdsState();
    explicit PdsState(PowerGrid&& graph);

    explicit PdsState(const PowerGrid& graph);

    PdsState(const PdsState&) = default;
    PdsState(PdsState&&) = default;
    PdsState& operator=(const PdsState&) = default;
    PdsState& operator=(PdsState&&) = default;

    void addEdge(Vertex source, Vertex target);
    void removeVertex(Vertex v);

    inline bool isZeroInjection(Vertex vertex) const { return m_graph[vertex].zero_injection; }

    inline PmuState activeState(Vertex vertex) const { return m_active.at(vertex); }

    inline bool isActive(Vertex vertex) const { return activeState(vertex) == PmuState::Active; }

    inline bool isBlank(Vertex vertex) const { return activeState(vertex) == PmuState::Blank; }

    inline bool isInactive(Vertex vertex) const { return activeState(vertex) == PmuState::Inactive; }

    inline bool isObserved(Vertex vertex) const { return m_observed.contains(vertex); }

    void createCheckpoint();

    std::span<Vertex> observedSinceCheckpoint();

    void restoreLastCheckpoint();

    inline size_t unobservedDegree(Vertex v) const {
        assert(m_unobserved_degree.at(v) == ranges::distance(m_graph.neighbors(v) | ranges::views::filter([this](auto v) { return !isObserved(v);})));
        return m_unobserved_degree.at(v);
    }

    inline const set<Vertex>& observed() const { return m_observed; }

    inline const map<Vertex, PmuState>& active() const { return m_active; }

    bool setActive(Vertex vertex);

    bool setInactive(Vertex vertex);

    bool setBlank(Vertex vertex);

    inline bool allObserved() const {
        return m_observed.size() == graph().numVertices();
    }

    size_t numObserved() const {
        return m_observed.size();
    }

    size_t numActive() const {
        return ranges::distance(graph().vertices() | ranges::views::filter([this](auto v) { return isActive(v); }));
    }

    size_t numInactive() const {
        return ranges::distance(graph().vertices() | ranges::views::filter([this](auto v) { return isInactive(v); }));
    }

    inline const PowerGrid& graph() const { return m_graph; }

    inline PowerGrid && graph() { return std::move(m_graph); }

    bool collapseLeaves();

    bool disableLowDegree();

    bool reduceObservedNonZi();

    bool collapseDegreeTwo();

    bool disableObservationNeighborhood();

    bool activateNecessaryNodes();

    bool collapseObservedEdges();

    bool solveTrivial();

    std::vector<PdsState> subproblems(bool nonZiSeparators = false) const;

};

void dominate(const PowerGrid& graph, const map<PowerGrid::vertex_descriptor, PmuState>& active, set<PowerGrid::vertex_descriptor>& observed);

bool propagate(const PowerGrid& graph, set<PowerGrid::vertex_descriptor>& observed, size_t max_unobserved);

bool observed(const PowerGrid& graph, const set<PowerGrid::vertex_descriptor> & observed);

} // namespace pds

#endif //PDS_PDS_HPP
