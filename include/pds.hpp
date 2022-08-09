//
// Created by max on 19.07.22.
//

#ifndef PDS_PDS_HPP
#define PDS_PDS_HPP

#include <range/v3/all.hpp>

#include <iostream>
#include <fstream>
#include <optional>

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

class PdsState {
public:
    using Vertex = PowerGrid::vertex_descriptor;
private:
    pds::map<Vertex, int> m_unobserved_degree;
    pds::set<Vertex> m_deleted;
    pds::set<Vertex> m_observed;
    pds::map<Vertex, PmuState> m_active;
    PowerGrid m_graph;

    inline void propagate(Vertex vertex) {
        if (is_observed(vertex) && zero_injection(vertex) && m_unobserved_degree[vertex] == 1) {
            for (auto w: m_graph.neighbors(vertex)) {
                observe(w);
            }
        }
    }

    void observe(Vertex vertex) {
        if (!is_observed(vertex)) {
            m_observed.insert(vertex);
            for (auto w: m_graph.neighbors(vertex)) {
                m_unobserved_degree[w] -= 1;
                assert(m_unobserved_degree[w] >= 0);
                propagate(w);
            }
            propagate(vertex);
        }
    }

    inline bool disableLowDegreeRecursive(
            PowerGrid::vertex_descriptor start,
            set<PowerGrid::vertex_descriptor>& seen
    ) {
        bool changed = false;
        seen.insert(start);
        if (m_graph.degree(start) <= 2 && zero_injection(start) && m_active.at(start) != PmuState::Active) {
            setInactive(start);
            changed = true;
        }
        for (auto w: m_graph.neighbors(start)) {
            if (!seen.contains(w)) {
                changed |= disableLowDegreeRecursive(w, seen);
            }
        }
        return changed;
    }

public:
    PdsState() = default;
    PdsState(PowerGrid&& graph) : m_graph(graph) {
        for (auto v: graph.vertices()) {
            m_unobserved_degree.emplace(v, m_graph.degree(v));
            m_active.emplace(v, PmuState::Blank);
        }
    }

    PdsState(const PowerGrid& graph) : m_graph(graph) {
        for (auto v: graph.vertices()) {
            m_unobserved_degree.emplace(v, m_graph.degree(v));
            m_active.emplace(v, PmuState::Blank);
        }
    }

    inline bool zero_injection(Vertex vertex) const {
        return m_graph[vertex].zero_injection;
    }

    inline PmuState active_state(Vertex vertex) const {
        return m_active.at(vertex);
    }

    inline bool is_observed(Vertex vertex) const {
        return m_observed.contains(vertex);
    }

    inline const set<Vertex>& observed() const {
        return m_observed;
    }

    inline const map<Vertex, PmuState> active() const {
        return m_active;
    }

    bool setActive(Vertex vertex) {
        if (m_active.at(vertex) != PmuState::Active) {
            m_active[vertex] = PmuState::Active;
            observe(vertex);
            for (auto w: m_graph.neighbors(vertex)) {
                observe(w);
            }
        }
        return m_observed.size() == m_graph.numVertices();
    }

    bool setInactive(Vertex vertex) {
        assert(active_state(vertex) != PmuState::Active);
        if (active_state(vertex) != PmuState::Inactive) {
            m_active[vertex] = PmuState::Inactive;
            assert(active_state(vertex) == PmuState::Inactive);
            return true;
        } else {
            return false;
        }
    }

    inline bool allObserved() const {
        return ranges::all_of(m_graph.vertices(), [this](auto v) { return m_observed.contains(v); });
    }

    inline const PowerGrid& graph() const {
        return m_graph;
    }

    inline PowerGrid && graph() {
        return std::move(m_graph);
    };

    inline bool collapseLeaves() {
        bool changed = false;
        auto vertices = m_graph.vertices() | ranges::to<std::vector>();
        for (auto v: vertices) {
            if (m_graph.degree(v) == 1
                && m_active.at(v) != PmuState::Active
            ) {
                Vertex neighbor;
                for (auto w: m_graph.neighbors(v)) {
                    neighbor = w;
                }
                if (m_graph.degree(neighbor) == 2 && m_graph[neighbor].zero_injection) {
                    if (!zero_injection(v)) {
                        m_graph[neighbor].zero_injection = false;
                    }
                } else {
                    if (m_graph[neighbor].zero_injection) {
                        m_graph[neighbor].zero_injection = false;
                    } else {
                        setActive(neighbor);
                    }
                }
                m_graph.removeVertex(v);
                changed = true;
            } else if (m_unobserved_degree.at(v) == 0 && m_active.at(v) == PmuState::Blank && is_observed(v)) {
                setInactive(v);
                changed = true;
            }
        }
        return changed;
    }

    inline bool disableLowDegree() {
        bool changed = false;
        set<PowerGrid::vertex_descriptor> seen;
        for (auto v: m_graph.vertices()) {
            if (m_graph.degree(v) >= 3) {
                changed = disableLowDegreeRecursive(v, seen);
            }
        }
        return changed;
    }

    inline bool collapseDegreeTwo() {
        bool changed = false;
        auto vertices = m_graph.vertices() | ranges::to<std::vector>();
        for (auto v: vertices) {
            if (m_graph.degree(v) == 2 && zero_injection(v)) {
                std::vector<Vertex> neighbors = m_graph.neighbors(v) | ranges::to<std::vector>();
                auto [x, y] = std::tie(neighbors[0], neighbors[1]);
                if (
                        !m_graph.edge(x, y)
                        && ((zero_injection(x) && m_graph.degree(x) <= 2)
                        || (zero_injection(y) && m_graph.degree(y) <= 2))) {
                    m_graph.addEdge(x, y);
                    m_graph.removeVertex(v);
                }
            }
        }
        return changed;
    }
};

inline auto dominate(const PowerGrid& graph, const map<PowerGrid::vertex_descriptor, PmuState>& active, set<PowerGrid::vertex_descriptor>& observed) {
    for (auto v: graph.vertices()) {
        if (active.at(v) == PmuState::Active) {
            observed.insert(v);
            for (auto w: graph.neighbors(v)) {
                observed.insert(w);
            }
        }
    }
    return observed;
}

inline bool propagate(const PowerGrid& graph, set<PowerGrid::vertex_descriptor>& observed, size_t max_unobserved = 1) {
    pds::map<PowerGrid::vertex_descriptor, size_t> unobserved_degree;
    std::vector<PowerGrid::vertex_descriptor> queue;
    for (const auto& v: graph.vertices()) {
        unobserved_degree[v] = 0;
        for (const auto& w: graph.neighbors(v)) {
            unobserved_degree[v] += !observed.contains(w);
        }
        if (graph[v].zero_injection && observed.contains(v) && unobserved_degree.at(v) > 0 && unobserved_degree.at(v) <= max_unobserved) {
            queue.push_back(v);
        }
    }
    while (!queue.empty()) {
        auto v = queue.back();
        queue.pop_back();
        for (const auto& w: graph.neighbors(v)) {
            if (!observed.contains(w)) {
                observed.insert(w);
                if (graph[w].zero_injection && observed.contains(w) && unobserved_degree.at(w) > 0 && unobserved_degree.at(w) <= max_unobserved) {
                    queue.push_back(w);
                }
                for (const auto& u: graph.neighbors(w)) {
                    unobserved_degree[u] -= 1;
                    auto deg = unobserved_degree.at(u);
                    if (graph[u].zero_injection && observed.contains(u) && deg > 0 && deg == max_unobserved) {
                        queue.push_back(u);
                    }
                }
            }
        }
    }
    return true;
}

inline bool observed(const PowerGrid& graph, const set<PowerGrid::vertex_descriptor> & observed) {
    return ranges::all_of(graph.vertices(), [&] (const auto& v) {
        return observed.contains(v);
    });
}

} // namespace pds

#endif //PDS_PDS_HPP
