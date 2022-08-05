//
// Created by max on 19.07.22.
//

#ifndef PDS_PDS_HPP
#define PDS_PDS_HPP

#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/graph/connected_components.hpp>

#include <range/v3/all.hpp>

#include <iostream>
#include <fstream>
#include <optional>

#include "utility.hpp"
#include "map.hpp"
#include "graph.hpp"

template<class T>
struct PrintType;

template<class T>
void print_type(T&&) { PrintType<T>{}; }

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

using PowerGrid = pds::AdjGraph<boost::listS, boost::setS, boost::undirectedS, Bus>;//boost::adjacency_list<boost::listS, boost::setS, boost::undirectedS, Bus>;

PowerGrid import_graphml(const std::string& filename, bool all_zero_injection = false);

class PdsState : public PowerGrid {
public:
    using Vertex = PowerGrid::vertex_descriptor;
    using Edge = PowerGrid::edge_descriptor;
    using VertexID = decltype(Bus::id);
private:
    pds::map<VertexID, int> m_unobserved_degree;
    pds::set<VertexID> m_deleted;
    pds::set<VertexID> m_observed;
    pds::map<VertexID, PmuState> m_active;

    void propagate(Vertex vertex) {
        if (is_observed(vertex) && zero_injection(vertex) && m_unobserved_degree[vertex_id(vertex)] == 1) {
            for (auto w: adjacent_vertices(vertex)) {
                observe(w);
            }
        }
    }

    void observe(Vertex vertex) {
        if (!is_observed(vertex)) {
            m_observed.insert(vertex_id(vertex));
            for (auto w: adjacent_vertices(vertex)) {
                m_unobserved_degree[vertex_id(w)] -= 1;
                assert(m_unobserved_degree[vertex_id(w)] >= 0);;
                propagate(w);
            }
        }
    }
public:
    PdsState() = default;
    PdsState(PowerGrid&& graph) : PowerGrid(graph) {
        for (auto v: graph.vertices()) {
            m_unobserved_degree.emplace(vertex_id(v), degree(v));
            m_active.emplace(vertex_id(v), PmuState::Blank);
        }
    }

    PdsState(const PowerGrid& graph) : PowerGrid(graph) {
        for (auto v: graph.vertices()) {
            m_unobserved_degree.emplace(vertex_id(v), degree(v));
            m_active.emplace(vertex_id(v), PmuState::Blank);
        }
    }

    inline VertexID vertex_id(Vertex vertex) const {
        return (*this)[vertex].id;
    }

    inline bool zero_injection(Vertex vertex) const {
        return (*this)[vertex].zero_injection;
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
            for (auto w: adjacent_vertices(vertex)) {
                observe(w);
            }
        }
        return m_observed.size() == num_vertices();
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
        return ranges::all_of(vertices(), [this](auto v) { return m_observed.contains(vertex_id(v)); });
    }

    inline const PowerGrid& graph() const {
        return *this;
    }
};

namespace detail {
template<class ZeroInjectionProperty, class ActivePropertyMap, class ObservedPropertyMap>
class Presolver {
public:
    using Vertex = PowerGrid::vertex_descriptor;
private:
    pds::map<Vertex, size_t> unobserved_degree_map;
    boost::associative_property_map<decltype(unobserved_degree_map)> unobserved_degree;
    pds::set<Vertex> deleted;
    pds::map<Vertex, Vertex> parent;
    PowerGrid& graph;
    ZeroInjectionProperty& zero_injection;
    ActivePropertyMap& active;
    ObservedPropertyMap& observed;


    void handle_degree_one(const Vertex& vertex) {
        assert(unobserved_degree[vertex] == 1);
        for (const auto& w: graph.adjacent_vertices(vertex)) {
            if (unobserved_degree[w] > 2) {
                if (zero_injection[w]) {
                    set_active(w);
                } else {
                    mark_deleted(vertex);
                    set_inactive(vertex);
                    zero_injection[w] = true;
                }
            }
        }
    }

    Vertex recurseDeg2Path(const Vertex& start, std::vector<Vertex>& path) {
        auto neigh = graph.adjacent_vertices(start)
                     | ranges::views::filter([this, &path](const auto& w) { return path.back() != w; })
                     | ranges::to<std::vector>();
        if (neigh.size() == 1 && zero_injection[start]) { // end of a path
            path.push_back(start);
            return recurseDeg2Path(neigh[0], path);
        } else if (neigh.size() == 0){
            path.push_back(start);
            return start;
        } else {
            return start;
        }
    }

    void handle_path(const Vertex& vertex) {
        if (unobserved_degree[vertex] <= 2) return;
        auto deg2neigh = graph.adjacent_vertices(vertex)
                         | ranges::views::filter([this](const auto& w) { return unobserved_degree[w] <= 2 && zero_injection[w]; })
                         | ranges::to<std::vector>();
        for (auto neigh: deg2neigh) {
            std::vector<Vertex> path = {vertex};
            auto last = recurseDeg2Path(neigh, path);
            std::cout << "found path of length " << path.size() << std::endl;
            if (last == vertex) {
                set_active(vertex);
                for (size_t i = 1; i < path.size(); ++i) {
                    set_inactive(path[i]);
                }
                path.pop_back();
            }
            for (size_t i = 1; i < path.size() - 1; ++i) {
                mark_deleted(path[i]);
                set_inactive(path[i]);
            }
            if (path.back() != vertex && !graph.edge(vertex, path.back())) {
                graph.template add_edge(vertex, path.back());
                if (!observed[vertex]) add_unobserved_neighbor(path.back());
                if (!observed[path.back()]) add_unobserved_neighbor(vertex);
            }
        }
    }

    inline void remove_unobserved_neighbor(const Vertex& vertex) {
        unobserved_degree[vertex] -= 1;
    }

    inline void add_unobserved_neighbor(const Vertex& vertex) {
        unobserved_degree[vertex] += 1;
    }

    inline void mark_deleted(const Vertex& vertex) {
        deleted.insert(vertex);
        observe(vertex);
        //for (const auto& u: range_pair(adjacent_vertices(vertex, graph))) {
        //    for (const auto& v: range_pair(adjacent_vertices(vertex, graph))) {
        //        if (u != v && !boost::edge(u, v, graph).second) {
        //            boost::add_edge(u, v, graph);
        //            if (!observed[u]) {
        //                add_unobserved_neighbor(v);
        //            }
        //            if (!observed[v]) {
        //                add_unobserved_neighbor(u);
        //            }
        //        }
        //    }
        //}
        //clear_vertex(vertex, graph);
    }

    void observe(const Vertex& vertex) {
        if (!observed[vertex]) {
            for (const auto &w: graph.adjacent_vertices(vertex)) {
                remove_unobserved_neighbor(w);
            }
            observed[vertex] = true;
        }
    }

    void set_inactive(const Vertex& vertex) {
        assert(active[vertex] != PmuState::Active);
        if (active[vertex] == PmuState::Blank) {
            active[vertex] = PmuState::Inactive;
        }
    }

    void set_active(const Vertex& vertex) {
        assert(active[vertex] != PmuState::Inactive);
        if (active[vertex] == PmuState::Blank) {
            active[vertex] = PmuState::Active;
            for (const auto& w: graph.adjacent_vertices(vertex)) {
                if (!observed[vertex]) remove_unobserved_neighbor(w);
                observe(w);
            }
            observed[vertex] = true;
        }
    }

public:
    Presolver(
            PowerGrid& graph,
            ZeroInjectionProperty& zero_injection,
            ActivePropertyMap& active,
            ObservedPropertyMap& observed
    ) : graph(graph), zero_injection(zero_injection), active(active), observed(observed),
        unobserved_degree_map(), unobserved_degree(unobserved_degree_map), deleted(),
        parent(ranges::to<decltype(parent)>(
                graph.vertices()
                | ranges::views::transform([](auto v){return std::make_pair(v, v);})
        ))
    { }

    void run() {
        pds::map<Vertex, size_t> component_map;
        pds::map<Vertex, unsigned char> color_map;
        size_t num_components = boost::connected_components(
                graph,
                boost::associative_property_map(component_map),
                boost::color_map(boost::associative_property_map(color_map)));
        assert(num_components == 1); // for now we don't support multiple components
        for (const auto& v: graph.vertices()) {
            unobserved_degree[v] = ranges::distance(
                    graph.adjacent_vertices(v) | ranges::views::filter([this](const auto& w){ return !deleted.contains(w) && !observed[w]; })
            );
        }
        for (const auto& v: range_pair(vertices(graph))) {
            if (unobserved_degree[v] == 1) {
                //handle_degree_one(v);
            }
            handle_path(v);
        }
    }
};
}

template<class ZeroInjectionProperty, class ActivePropertyMap, class ObservedPropertyMap>
void presolve(PowerGrid& graph, ZeroInjectionProperty& zero_injection, ActivePropertyMap& active, ObservedPropertyMap& observed) {
    detail::Presolver presolver(graph, zero_injection, active, observed);
    presolver.run();
}

template<typename T> struct PrintType;
template<class Graph, class ActivePropertyMap, class ObservedPropertyMap>
ObservedPropertyMap dominate(const Graph& graph, const ActivePropertyMap& active, ObservedPropertyMap observed) {
    for (auto v: graph.vertices()) {
        if (active[v] == PmuState::Active) {
            observed[v] = true;
            auto adjacent = boost::adjacent_vertices(v, graph);
            for (auto it = adjacent.first; it != adjacent.second; ++it) {
                observed[*it] = true;
            }
        }
    }
    return observed;
}

template<class ObservedPropertyMap>
std::vector<PowerGrid::vertex_descriptor> propagation_targets_k(
        PowerGrid& graph,
        const PowerGrid::vertex_descriptor& vertex,
        ObservedPropertyMap& observed,
        size_t max_unobserved = 1)
{
    if (!observed[vertex]) return {};
    std::vector<PowerGrid::vertex_descriptor> targets;
    for (const auto& w: graph.adjacent_vertices(vertex)) {
        if (!observed[w]) {
            if (targets.size() >= max_unobserved) {
                return {};
            } else {
                targets.push_back(w);
            }
        }
    }
    return targets;
}

template<class ActivePropertyMap>
std::vector<PowerGrid::vertex_descriptor> activeVertices(const PowerGrid& graph, const ActivePropertyMap& active) {
    return ranges::to<std::vector>(graph.vertices() | ranges::views::filter([&active](const auto& v) { return active[v] == pds::PmuState::Active; }));
}

template<class ZeroInjectionProperty, class ObservedPropertyMap>
bool propagate(const PowerGrid& graph, const ZeroInjectionProperty& zero_injection, ObservedPropertyMap observed, size_t max_unobserved = 1) {
    pds::map<PowerGrid::vertex_descriptor, size_t> unobserved_degree_map;
    auto unobserved_degree = boost::associative_property_map(unobserved_degree_map);
    std::vector<PowerGrid::vertex_descriptor> queue;
    for (const auto& v: graph.vertices()) {
        for (const auto& w: graph.adjacent_vertices(v)) {
            unobserved_degree[v] += !observed[w];
        }
        if (zero_injection[v] && observed[v] && unobserved_degree[v] > 0 && unobserved_degree[v] <= max_unobserved) {
            queue.push_back(v);
        }
    }
    while (!queue.empty()) {
        auto v = queue.back();
        queue.pop_back();
        for (const auto& w: graph.adjacent_vertices(v)) {
            if (!observed[w]) {
                observed[w] = true;
                if (zero_injection[w] && observed[w] && unobserved_degree[w] > 0 && unobserved_degree[w] <= max_unobserved) {
                    queue.push_back(w);
                }
                for (const auto& u: graph.adjacent_vertices(w)) {
                    unobserved_degree[u] -= 1;
                    auto deg = unobserved_degree[u];
                    if (zero_injection[u] && observed[u] && deg > 0 && deg == max_unobserved) {
                        queue.push_back(u);
                    }
                }
            }
        }
    }
    return true;
}

template<class ObservedPropertyMap>
bool observed(const PowerGrid& graph, const ObservedPropertyMap& observed) {
    return ranges::all_of(graph.vertices(), [&] (const auto& v) {
        return observed[v];
    });
}

} // namespace pds

#endif //PDS_PDS_HPP
