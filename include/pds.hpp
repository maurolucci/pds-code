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
#include <setgraph/graph.hpp>

namespace pds {

enum class PmuState {
    Blank,
    Active,
    Inactive
};

struct Bus {
    std::string name;
    long id;
    bool zero_injection;
};

using PowerGrid = setgraph::SetGraph<Bus, setgraph::Empty, setgraph::EdgeDirection::Undirected>;

PowerGrid import_graphml(const std::string& filename, bool all_zero_injection = false);

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
    pds::map<Vertex, int> m_unobserved_degree;
    pds::set<Vertex> m_deleted;
    pds::set<Vertex> m_observed;
    pds::map<Vertex, PmuState> m_active;
    PowerGrid m_graph;

    void propagate(Vertex vertex);

    void observe(Vertex vertex);

    inline bool disableLowDegreeRecursive(
            PowerGrid::vertex_descriptor start,
            set<PowerGrid::vertex_descriptor>& seen
    );

public:
    PdsState() = default;
    explicit PdsState(PowerGrid&& graph) : m_graph(graph) {
        for (auto v: graph.vertices()) {
            m_unobserved_degree.emplace(v, m_graph.degree(v));
            m_active.emplace(v, PmuState::Blank);
        }
    }

    explicit PdsState(const PowerGrid& graph) : m_graph(graph) {
        for (auto v: graph.vertices()) {
            m_unobserved_degree.emplace(v, m_graph.degree(v));
            m_active.emplace(v, PmuState::Blank);
        }
    }

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

    inline const set<Vertex>& observed() const { return m_observed; }

    inline const map<Vertex, PmuState>& active() const { return m_active; }

    bool setActive(Vertex vertex);

    bool setInactive(Vertex vertex);

    bool setBlank(Vertex vertex);

    inline bool allObserved() const {
        return ranges::all_of(m_graph.vertices(), [this](auto v) { return m_observed.contains(v); });
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
