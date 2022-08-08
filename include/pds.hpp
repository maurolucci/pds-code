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
    using VertexID = decltype(Bus::id);
private:
    pds::map<VertexID, int> m_unobserved_degree;
    pds::set<VertexID> m_deleted;
    pds::set<VertexID> m_observed;
    pds::map<VertexID, PmuState> m_active;
    PowerGrid m_graph;

    inline void propagate(Vertex vertex) {
        if (is_observed(vertex) && zero_injection(vertex) && m_unobserved_degree[vertex_id(vertex)] == 1) {
            for (auto w: m_graph.neighbors(vertex)) {
                observe(w);
            }
        }
    }

    void observe(Vertex vertex) {
        if (!is_observed(vertex)) {
            m_observed.insert(vertex_id(vertex));
            for (auto w: m_graph.neighbors(vertex)) {
                m_unobserved_degree[vertex_id(w)] -= 1;
                assert(m_unobserved_degree[vertex_id(w)] >= 0);;
                propagate(w);
            }
            propagate(vertex);
        }
    }
public:
    PdsState() = default;
    PdsState(PowerGrid&& graph) : m_graph(graph) {
        for (auto v: graph.vertices()) {
            m_unobserved_degree.emplace(vertex_id(v), m_graph.degree(v));
            m_active.emplace(vertex_id(v), PmuState::Blank);
        }
    }

    PdsState(const PowerGrid& graph) : m_graph(graph) {
        for (auto v: graph.vertices()) {
            m_unobserved_degree.emplace(vertex_id(v), m_graph.degree(v));
            m_active.emplace(vertex_id(v), PmuState::Blank);
        }
    }

    inline VertexID vertex_id(Vertex vertex) const {
        return m_graph[vertex].id;
    }

    inline bool zero_injection(Vertex vertex) const {
        return m_graph[vertex].zero_injection;
    }

    inline PmuState active_state(Vertex vertex) const {
        return m_active.at(vertex_id(vertex));
    }

    inline bool is_observed(Vertex vertex) const {
        return m_observed.contains(vertex_id(vertex));
    }

    bool set_active(Vertex vertex) {
        if (m_active.at(vertex_id(vertex)) != PmuState::Active) {
            m_active.emplace(vertex_id(vertex), PmuState::Active);
            observe(vertex);
            for (auto w: m_graph.neighbors(vertex)) {
                observe(w);
            }
        }
        return m_observed.size() == m_graph.numVertices();
    }

    bool set_inactive(Vertex vertex) {
        assert(active_state(vertex) != PmuState::Active);
        if (active_state(vertex) != PmuState::Inactive) {
            m_active.emplace(vertex_id(vertex), PmuState::Inactive);
            return true;
        } else {
            return false;
        }
    }

    inline bool all_observed() const {
        return ranges::all_of(m_graph.vertices(), [this](auto v) { return m_observed.contains(vertex_id(v)); });
    }

    inline const PowerGrid& graph() const {
        return m_graph;
    }

    inline PowerGrid && graph() {
        return std::move(m_graph);
    };
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
