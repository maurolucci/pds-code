//
// Created by max on 01.08.22.
//

#include <setgraph/boost_adapter.hpp>
#include "pds.hpp"
#include "graphml/graphml.hpp"

#include <boost/property_map/function_property_map.hpp>

namespace pds {

PowerGrid import_graphml(const std::string& filename, bool all_zero_injection) {
    PowerGrid graph;
    boost::dynamic_properties attr(boost::ignore_other_properties);
    map<PowerGrid::vertex_descriptor, long> zero_injection_data;
    boost::associative_property_map zero_injection(zero_injection_data);
    auto id_map = [&graph](const PowerGrid::vertex_descriptor& vertex) -> long& { return graph[vertex].id;};
    auto name_map = [&graph](PowerGrid::vertex_descriptor vertex) -> std::string& { return graph[vertex].name; };
    boost::function_property_map<decltype(id_map), PowerGrid::vertex_descriptor> id(id_map);
    boost::function_property_map<decltype(name_map), PowerGrid::vertex_descriptor> name(name_map);
    attr.property("zero_injection", zero_injection);//make_vector_property_map(long_zi, graph));
    attr.property("name", name);
    attr.property("id", id);
    std::ifstream graph_in(filename);
    pds::read_graphml(graph_in, graph, attr);
    graph_in.close();
    for (auto v: graph.vertices()) {
        graph[v].zero_injection = zero_injection[v] || all_zero_injection;
    }
    return graph;
}

set<PowerGrid::vertex_descriptor> observationNeighborhood(const PowerGrid &graph, const set<PowerGrid::vertex_descriptor> &starts) {
    set<PowerGrid::vertex_descriptor> observed;
    for (auto v: starts) {
        observed.insert(v);
        for (auto w: graph.neighbors(v)) {
            observed.insert(w);
        }
    }
    propagate(graph, observed);
    return observed;
}

void PdsState::addEdge(Vertex source, Vertex target) {
    if (!m_graph.edge(source, target)) {
        m_graph.addEdge(source, target);
        if (!isObserved(source)) {
            m_unobserved_degree[target] += 1;
        }
        if (!isObserved(target)) {
            m_unobserved_degree[source] += 1;
        }
    }
}

void PdsState::removeVertex(Vertex v) {
    if (!isObserved(v)) {
        for (auto w: m_graph.neighbors(v)) {
            m_unobserved_degree[w] -= 1;
            propagate(w);
        }
    }
    m_graph.removeVertex(v);
}

void PdsState::propagate(PdsState::Vertex vertex) {
    if (isObserved(vertex) && isZeroInjection(vertex) && m_unobserved_degree[vertex] == 1) {
        for (auto w: m_graph.neighbors(vertex)) {
            observe(w);
        }
    }
}

void PdsState::observe(PdsState::Vertex vertex) {
    if (!isObserved(vertex)) {
        m_observed.insert(vertex);
        for (auto w: m_graph.neighbors(vertex)) {
            m_unobserved_degree[w] -= 1;
            assert(m_unobserved_degree[w] >= 0);
            propagate(w);
        }
        propagate(vertex);
    }
}

bool PdsState::disableLowDegreeRecursive(PdsState::Vertex start, set<PdsState::Vertex> &seen) {
    auto hasBlankNeighbor = [this] (auto v) {
        return ranges::any_of(m_graph.neighbors(v), [this](auto w) { return isActive(w) || isBlank(w);} );
    };
    bool changed = false;
    seen.insert(start);
    if (m_graph.degree(start) <= 2 && hasBlankNeighbor(start) && isZeroInjection(start) && isBlank(start)) {
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

bool PdsState::setBlank(PdsState::Vertex vertex) {
    assert(!isActive(vertex));
    m_active[vertex] = PmuState::Blank;
    return m_observed.size() == m_graph.numVertices();
}

bool PdsState::setActive(PdsState::Vertex vertex) {
    m_active[vertex] = PmuState::Active;
    observe(vertex);
    for (auto w: m_graph.neighbors(vertex)) {
        observe(w);
    }
    return m_observed.size() == m_graph.numVertices();
}

bool PdsState::setInactive(PdsState::Vertex vertex) {
    assert(!isActive(vertex));
    if (!isInactive(vertex)) {
        m_active[vertex] = PmuState::Inactive;
        assert(activeState(vertex) == PmuState::Inactive);
        return true;
    } else {
        return false;
    }
}

bool PdsState::collapseLeaves() {
    bool changed = false;
    auto vertices = m_graph.vertices() | ranges::to<std::vector>();
    for (auto v: vertices) {
        if (m_graph.degree(v) == 1 && !isActive(v)) {
            Vertex neighbor;
            for (auto w: m_graph.neighbors(v)) {
                neighbor = w;
            }
            if (m_graph.degree(neighbor) == 2 && isZeroInjection(neighbor)) {
                if (!isZeroInjection(v)) {
                    m_graph[neighbor].zero_injection = false;
                }
            } else {
                if (isZeroInjection(neighbor)) {
                    m_graph[neighbor].zero_injection = false;
                } else {
                    setActive(neighbor);
                }
            }
            removeVertex(v);
            changed = true;
        } else if (m_unobserved_degree.at(v) == 0 && isBlank(v) && isObserved(v)) {
            setInactive(v);
            changed = true;
        } else if (m_graph.degree(v) == 0 && !isObserved(v) && isBlank(v)) {
            setActive(v);
        }
    }
    return changed;
}

bool PdsState::disableLowDegree() {
    bool changed = false;
    set<PdsState::Vertex> seen;
    for (auto v: m_graph.vertices()) {
        if (m_graph.degree(v) >= 3) {
            changed = disableLowDegreeRecursive(v, seen);
        }
    }
    return changed;
}

bool PdsState::reduceObservedNonZi() {
    bool changed = false;
    auto vertices = m_graph.vertices() | ranges::to<std::vector>();
    for (auto v: vertices) {
        if (isObserved(v) && isInactive(v) && !isZeroInjection(v)) {
            removeVertex(v);
            changed = true;
        }
    }
    return changed;
}

bool PdsState::collapseDegreeTwo() {
    bool changed = false;
    auto vertices = m_graph.vertices() | ranges::to<std::vector>();
    for (auto v: vertices) {
        if (m_graph.degree(v) == 2 && isZeroInjection(v)) {
            std::vector<Vertex> neighbors = m_graph.neighbors(v) | ranges::to<std::vector>();
            auto [x, y] = std::tie(neighbors[0], neighbors[1]); { }
            if (!m_graph.edge(x, y)) {
                if (isInactive(v)) {
                    if((isZeroInjection(x) && m_graph.degree(x) <= 2) || (isZeroInjection(y) && m_graph.degree(y) <= 2)) {
                        addEdge(x, y);
                        removeVertex(v);
                        changed = true;
                    }
                }
                for (auto [s, t]: {std::tie(x, y), std::tie(y, x)}) {
                    if (!isZeroInjection(s) && isInactive(s) && !isZeroInjection(t) && isBlank(t)) {
                        setActive(t);
                        changed = true;
                    }
                }
            } else {
                if (m_graph.degree(x) == 2 && isBlank(y)) {
                    setActive(y);
                    changed = true;
                } else if (m_graph.degree(y) == 2 && isBlank(x)) {
                    setActive(x);
                    changed = true;
                }
            }
        }
    }
    return changed;
}

bool PdsState::disableObservationNeighborhood() {
    bool changed = false;
    map<Vertex, set<Vertex>> cachedNeighborhood;
    set<Vertex> active;
    for (auto v: m_graph.vertices()) {
        if (isActive(v)) active.insert(v);
    }
    auto activeNeighborhood = observationNeighborhood(m_graph, active);
    for (auto v: m_graph.vertices()) {
        if (isBlank(v)) {
            if (isSuperset(activeNeighborhood, m_graph.neighbors(v))) {
                setInactive(v);
                changed = true;
            }
        }
    }
    auto vertices = m_graph.vertices() | ranges::to<std::vector>();
    ranges::sort(vertices, [this](auto left, auto right) -> bool { return m_graph.degree(left) > m_graph.degree(right);});
    for (auto v: vertices) {
        if (isBlank(v)) {
            active.insert(v);
            auto observation = observationNeighborhood(m_graph, active);
            active.erase(v);
            for (auto w: observation) {
                if (v != w && isBlank(w)) {
                    if (isSuperset(observation, m_graph.neighbors(w))) {
                        setInactive(w);
                        changed = true;
                    }
                }
            }
        }
    }
    return changed;
}

bool PdsState::activateNecessaryNodes() {
    bool changed = false;
    set<Vertex> blankActive;
    for (auto v: m_graph.vertices()) {
        if (!isInactive(v)) blankActive.insert(v);
    }
    for (auto v: m_graph.vertices()) {
        if (isBlank(v)) {
            blankActive.erase(v);
            auto observed = observationNeighborhood(m_graph, blankActive);
            if (ranges::any_of(m_graph.vertices(), [&observed](auto v) -> bool { return !observed.contains(v);})) {
                setActive(v);
                changed = true;
            }
            blankActive.insert(v);
        }
    }
    return changed;
}

map<PdsState::Vertex, PdsState::Vertex> findClosestActive(const PdsState& state) {
    using Vertex = PdsState::Vertex;
    std::deque<Vertex> queue;
    map<Vertex, Vertex> closestActive;
    for (auto v: state.graph().vertices()) {
        if (state.isActive(v)) {
            closestActive.emplace(v, v);
            queue.push_back(v);
        }
    }
    while (!queue.empty()) {
        auto v = queue.front();
        queue.pop_front();
        for (auto w: state.graph().neighbors(v)) {
            if (!closestActive.contains(w)) {
                closestActive.emplace(w, closestActive[v]);
                queue.push_back(w);
            }
        }
    }
    return closestActive;
}

bool PdsState::collapseObservedEdges() {
    bool changed = false;
    auto isObservedEdge = [this](PowerGrid::edge_descriptor e) {
        auto [s, t] = m_graph.endpoints(e);
        return isObserved(s) && isObserved(t);
    };
    auto edges = m_graph.edges() | ranges::views::filter(isObservedEdge) | ranges::to<std::vector>();
    if (edges.empty()) return false;
    auto closestActive = findClosestActive(*this);

    for (auto e: edges) {
        auto [s, t] = m_graph.endpoints(e);
        if (s == closestActive.at(t) || t == closestActive.at(s)) continue;
        m_graph.removeEdge(e);
        changed = true;
        auto hasObservedNeighbor = [this](Vertex v) -> bool { return ranges::any_of(m_graph.neighbors(v), [this](auto w) { return isActive(w);});};
        if (!isActive(s)) {
            if (!hasObservedNeighbor(s)) addEdge(s, closestActive.at(s));
        }
        if (!isActive(t)) {
            if (!hasObservedNeighbor(t)) addEdge(t, closestActive.at(t));
        }
    }
    return changed;
}

bool isBlockingVertex(const PdsState& state, PdsState::Vertex vertex) {
    return !state.isZeroInjection(vertex) && state.isObserved(vertex) && state.isInactive(vertex);
}

inline bool isBlockedEdge(const PdsState& state, PdsState::Vertex source, PdsState::Vertex target) {
    return state.isObserved(source) && state.isObserved(target)
        && !state.isActive(source) && !state.isActive(target)
        && !(state.isZeroInjection(source) && state.isZeroInjection(target));
}

void recurseComponent(
        PdsState::Vertex start,
        const PdsState& state,
        set<PdsState::Vertex>& seen,
        set<PdsState::Vertex>& currentComponent,
        bool nonZiSeparators
) {
    if (seen.contains(start)) return;
    seen.insert(start);
    currentComponent.insert(start);
    for (auto w: state.graph().neighbors(start)) {
        if (!state.isActive(w)) {
            if (!seen.contains(w) && (!nonZiSeparators || !isBlockedEdge(state, start, w))) {
                recurseComponent(w, state, seen, currentComponent, nonZiSeparators);
            }
        } else {
            currentComponent.insert(w);
        }
    }
}

bool PdsState::solveTrivial() {
    if (!allObserved()) {
        auto blank_filter = [this] (Vertex v) -> bool { return isBlank(v); };
        auto possibleLocations = m_graph.vertices() | ranges::views::filter(blank_filter) | ranges::to<std::vector>();
        if (possibleLocations.size() == 1) {
            setActive(possibleLocations[0]);
        }
    }
    return allObserved();
}

std::vector<PdsState> PdsState::subproblems(bool nonZiSeparators) const {
    std::vector<PdsState> components;
    set<Vertex> seen;
    for (auto v: m_graph.vertices()) {
        if (!isActive(v) && !seen.contains(v)) {
            set<Vertex> currentComponent;
            recurseComponent(v, *this, seen, currentComponent, nonZiSeparators);
            PowerGrid subgraph{graph()};
            for (auto x: subgraph.vertices()) {
                if (!currentComponent.contains(x)) {
                    subgraph.removeVertex(x);
                }
            }
            auto dummy = subgraph.addVertex();
            auto oneActive = dummy;
            PdsState subproblem(std::move(subgraph));
            for (auto x: subproblem.graph().vertices()) {
                if (x == dummy) continue;
                if (isActive(x)) {
                    subproblem.setActive(x);
                    oneActive = x;
                } else if (isInactive(x)) {
                    subproblem.setInactive(x);
                }
            }
            bool dummyNeeded = false;
            for (auto x: subproblem.graph().vertices()) {
                if (x == dummy) continue;
                if (isObserved(x) && !subproblem.isObserved(x)) {
                    subproblem.addEdge(x, oneActive);
                    dummyNeeded = true;
                }
            }
            subproblem.setActive(oneActive);
            if (oneActive != dummy || !dummyNeeded) {
                subproblem.removeVertex(dummy);
            }
            components.emplace_back(std::move(subproblem));
        }
    }
    return components;
}

void dominate(const PowerGrid &graph, const map<PowerGrid::vertex_descriptor, PmuState> &active, set<PowerGrid::vertex_descriptor> &observed) {
    for (auto v: graph.vertices()) {
        if (active.at(v) == PmuState::Active) {
            observed.insert(v);
            for (auto w: graph.neighbors(v)) {
                observed.insert(w);
            }
        }
    }
}

bool propagate(const PowerGrid &graph, set<PowerGrid::vertex_descriptor> &observed, size_t max_unobserved) {
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

bool observed(const PowerGrid &graph, const set<PowerGrid::vertex_descriptor> &observed) {
    return ranges::all_of(graph.vertices(), [&] (const auto& v) {
        return observed.contains(v);
    });
}

}
