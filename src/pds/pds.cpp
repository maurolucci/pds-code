//
// Created by max on 01.08.22.
//

#include "mpgraphs/boost_adapter.hpp"
#include "pds.hpp"
#include "graphml/graphml.hpp"

#include <boost/property_map/function_property_map.hpp>
#include "utility.hpp"

namespace pds {

//void validateState(const PdsState& state) {
//    auto active = state.active()
//            | ranges::views::filter([](const auto& xy) { return xy.second == PmuState::Active; })
//            | ranges::views::transform([](const auto& xy) { return xy.first; })
//            | ranges::to<set<PdsState::Vertex>>();
//    auto observed = observationNeighborhood(state.graph(), active);
//    for (auto v: state.graph().vertices()) {
//        unused(v);
//        assert(observed.contains(v) == state.isObserved(v));
//    }
//}

void exportGraphml(const PowerGrid& graph, std::ostream& out) {
    map<PowerGrid::VertexDescriptor, long> zero_injection_data;
    for (auto v: graph.vertices()) {
        if (graph[v].zero_injection) {
            zero_injection_data[v] = 1;
        }
    }
    boost::dynamic_properties attr(boost::ignore_other_properties);
    boost::associative_property_map zero_injection(zero_injection_data);
    auto id_map = [&graph](const PowerGrid::VertexDescriptor& vertex) -> long { return graph[vertex].id;};
    auto name_map = [&graph](PowerGrid::VertexDescriptor vertex) -> std::string { return graph[vertex].name; };
    boost::function_property_map<decltype(id_map), PowerGrid::VertexDescriptor> id(id_map);
    boost::function_property_map<decltype(name_map), PowerGrid::VertexDescriptor> name(name_map);
    attr.property("zero_injection", zero_injection);
    attr.property("name", name);
    attr.property("id", id);
    pds::write_graphml(out, graph, id, attr, true);
}

PdsState::PdsState() : PdsState(PowerGrid{}) { }

PdsState::PdsState(const pds::PowerGrid &graph) : PdsState(PowerGrid{graph}) { }

PdsState::PdsState(PowerGrid&& graph) : m_numActive{0}, m_numInactive{0}, m_graph(graph) {
    for (auto v: m_graph.vertices()) {
        m_graph.removeEdge(v, v);
        m_unobserved_degree[v] = m_graph.degree(v);
    }
    for (auto v: m_graph.vertices()) {
        switch (m_graph[v].pmu) {
            case PmuState::Active:
                m_graph[v].pmu = PmuState::Blank;
                setActive(v);
                break;
            case PmuState::Inactive:
                setInactive(v);
                break;
            default:
                break;
        }
    }
}

void PdsState::addEdge(Vertex source, Vertex target) {
    assert(source != target);
    if (!m_graph.edge(source, target)) {
        m_graph.addEdge(source, target);
        if (!isObserved(source)) {
            m_unobserved_degree[target] += 1;
        }
        if (!isObserved(target)) {
            m_unobserved_degree[source] += 1;
        }
    }
#ifndef NDEBUG
    for (auto v: {source, target}) {
        assert(m_unobserved_degree[v] == ranges::distance(
                m_graph.neighbors(v) | ranges::views::filter([this](auto v) { return !isObserved(v); })));
    }
#endif
}

void PdsState::removeVertex(Vertex v) {
    if (!isObserved(v)) {
        for (auto w: m_graph.neighbors(v)) {
            m_unobserved_degree[w] -= 1;
        }
    }
    assert(m_unobserved_degree[v] == ranges::distance(m_graph.neighbors(v) | ranges::views::filter([this](auto v) { return !isObserved(v); })));
    if (isActive(v)) --m_numActive;
    if (isInactive(v)) --m_numInactive;
    m_dependencies.removeVertex(v);
    m_unobserved_degree.erase(v);
    m_graph.removeVertex(v);
}

void PdsState::propagate(std::vector<Vertex> &queue) {
    while (!queue.empty()) {
        auto v = queue.back();
        queue.pop_back();
        if (isObserved(v) && isZeroInjection(v) && m_unobserved_degree[v] == 1) {
            for (auto w: m_graph.neighbors(v)) {
                if (!isObserved(w)) {
                    observeOne(w, v, queue);
                }
            }
        }
        assert(isActive(v) || !m_dependencies.hasVertex(v) || m_dependencies.inDegree(v) > 0);
    }
}

bool PdsState::observeOne(Vertex vertex, Vertex origin, std::vector<Vertex>& queue) {
    if (!isObserved(vertex)) {
        m_dependencies.getOrAddVertex(vertex);
        if (origin != vertex) m_dependencies.addEdge(origin, vertex);
        if (m_unobserved_degree[vertex] == 1) queue.push_back(vertex);
        for (auto w: m_graph.neighbors(vertex)) {
            m_unobserved_degree[w] -= 1;
            if (m_unobserved_degree[w] == 1 && isObserved(w) && isZeroInjection(w)) queue.push_back(w);
        }
        assert(isActive(vertex) || !m_dependencies.hasVertex(vertex) || m_dependencies.inDegree(vertex) > 0);
        return true;
    } else {
        return false;
    }
}

bool PdsState::observe(Vertex vertex, Vertex origin) {
    assert(isObserved(origin) || isActive(origin));
    assert(isActive(origin) || isZeroInjection(origin));
    std::vector<Vertex> queue;
    if (observeOne(vertex, origin, queue)) {
        propagate(queue);
        return true;
    } else {
        return false;
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
    if (isInactive(vertex)) {
        --m_numInactive;
        m_graph[vertex].pmu = PmuState::Blank;
    }
    return allObserved();
}

bool PdsState::setActive(PdsState::Vertex vertex) {
    assert(!isInactive(vertex));
    if (!isActive(vertex)) {
        ++m_numActive;
        m_graph[vertex].pmu = PmuState::Active;
        if (m_dependencies.hasVertex(vertex)) {
            while (m_dependencies.inDegree(vertex) > 0) {
                auto edge = *m_dependencies.inEdges(vertex).begin();
                m_dependencies.removeEdge(edge);
            }
        }
        std::vector<Vertex> queue;
        observeOne(vertex, vertex, queue);
        for (auto w: m_graph.neighbors(vertex)) {
            observeOne(w, vertex, queue);
        }
        propagate(queue);
    }
    return allObserved();
}

bool PdsState::unsetActive(PdsState::Vertex vertex) {
    assert(isActive(vertex));
    if (isActive(vertex)) {
        std::vector<Vertex> propagating;
        std::vector<Vertex> queue;
        set<Vertex> enqueued;

        m_graph[vertex].pmu = PmuState::Blank;
        --m_numActive;
        queue.push_back(vertex);
        enqueued.insert(vertex);
        while (!queue.empty()) {
            auto v = queue.back();
            queue.pop_back();
            assert(!isActive(v));
            std::optional<Vertex> observer;
            for (auto w: m_graph.neighbors(v)) {
                assert(m_unobserved_degree[w] == ranges::distance(
                        m_graph.neighbors(w) | ranges::views::filter([this](auto v) { return !isObserved(v); })));
                assert(isObserved(v));
                m_unobserved_degree[w] += 1;
                if (!isActive(w)) {
                    // w is observed from v
                    if (m_dependencies.edge(v, w).has_value()) {
                        if (!enqueued.contains(w)) {
                            queue.push_back(w);
                            enqueued.insert(w);
                        }
                    } else if (isObserved(w)) {
                        for (auto x: m_dependencies.neighbors(w)) {
                            if (!enqueued.contains(x)) {
                                queue.push_back(x);
                                enqueued.insert(x);
                            }
                        }
                        propagating.push_back(w);
                    }
                } else {
                    observer = {w};
                    assert(isActive(w));
                    assert(isObserved(w));
                }
            }
            // mark unobserved
            m_dependencies.removeVertex(v);
            if (observer.has_value()) {
                assert(isActive(*observer));
                assert(isObserved(*observer));
                observeOne(v, observer.value(), propagating);
            }
        }

        propagate(propagating);
    }
    return isObserved(vertex);
}

bool PdsState::setInactive(PdsState::Vertex vertex) {
    assert(!isActive(vertex));
    if (!isInactive(vertex)) {
        ++m_numInactive;
        m_graph[vertex].pmu = PmuState::Inactive;
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
            if (isInactive(v) || isBlank(neighbor)) {
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
            }
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
                    auto otherNeighbor = [this] (auto i, auto j) {
                        return (m_graph.neighbors(i) | ranges::views::filter([this,j](auto v) { return v != j && !isObserved(v);}) | ranges::to_vector)[0];
                    };
                    if (!isObserved(v) && isObserved(x) && isInactive(x) && isZeroInjection(x) && m_unobserved_degree[x] == 2 && m_graph.degree(otherNeighbor(x, v)) == 2 && isZeroInjection(otherNeighbor(x, v))) {
                        // v must be unobserved if x is observed
                        auto z = otherNeighbor(x, v);
                        if (!m_graph.edge(y, z).has_value()) {
                            assert(!isObserved(z) && z != v && m_graph.edge(v, x) && m_graph.edge(x, z) &&
                                   m_graph.degree(z) == 2);
                            assert(isObserved(x) && !isObserved(v) && !isObserved(z));
                            m_graph.removeEdge(x, z);
                            m_unobserved_degree[x] -= 1;
                            removeVertex(v);
                            addEdge(y, z);
                            changed = true;
                            assert(m_unobserved_degree[x] == 0);
                            assert(ranges::distance(m_graph.neighbors(x) |
                                                    ranges::views::filter([this](auto v) { return !isObserved(v); })) ==
                                   m_unobserved_degree[x]);
                        }
                    } else if (!isObserved(v) && isObserved(y) && isInactive(y) && isZeroInjection(y) && m_unobserved_degree[y] == 2 && m_graph.degree(otherNeighbor(y, v)) == 2 && isZeroInjection(otherNeighbor(y, v))) {
                        // v must be unobserved if x is observed
                        auto z = otherNeighbor(y, v);
                        if (!m_graph.edge(x, z).has_value()) {
                            assert(!isObserved(z) && z != v && m_graph.edge(v, y) && m_graph.edge(y, z) &&  m_graph.degree(z) == 2);
                            assert(isObserved(y) && !isObserved(v) && !isObserved(z));
                            m_graph.removeEdge(y, z);
                            m_unobserved_degree[y] -= 1;
                            removeVertex(v);
                            addEdge(x, z);
                            changed = true;
                            assert(m_unobserved_degree[y] == 0);
                            assert(ranges::distance(m_graph.neighbors(y) | ranges::views::filter([this](auto v) { return !isObserved(v);})) == m_unobserved_degree[y]);
                        }
                    } else if((isZeroInjection(x) && m_graph.degree(x) <= 2) || (isZeroInjection(y) && m_graph.degree(y) <= 2)) {
                        if (m_dependencies.edge(v, y).has_value()) {
                            m_dependencies.addEdge(x, y);
                        } else if (m_dependencies.edge(v, x).has_value()) {
                            m_dependencies.addEdge(y, x);
                        }
                        removeVertex(v);
                        addEdge(x, y);
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
    auto fullyObserved = [this](auto v) {
        return isObserved(v) && unobservedDegree(v) == 0;
    };
    std::vector<Vertex> blankVertices;
    std::vector<char> blankInactive;
    for (auto v: m_graph.vertices()) {
        if (isBlank(v)) {
            if (fullyObserved(v)) {
                setInactive(v);
                changed = true;
            } else {
                blankVertices.push_back(v);
                blankInactive.push_back(false);
            }
        }
    }
    ranges::sort(blankVertices, [this](auto left, auto right) -> bool { return m_graph.degree(left) > m_graph.degree(right);});
    size_t skipStart = 0;
    for (size_t i = 0; i < blankVertices.size(); ++i) {
        auto& v = blankVertices[i];
        auto& inactive = blankInactive[i];
        if (!inactive) {
            setActive(v);
            for (size_t j = skipStart; j < blankVertices.size(); ++j) {
                if (i != j) {
                    auto w = blankVertices[j];
                    if (!blankInactive[j] && fullyObserved(w)) {
                        setInactive(w);
                        changed = true;
                        if (j > i) {
                            std::swap(blankVertices[j], blankVertices.back());
                            std::swap(blankInactive[j], blankInactive.back());
                            blankVertices.pop_back();
                            blankInactive.pop_back();
                            --j; // j > 0 is guaranteed
                        } else {
                            std::swap(blankVertices[j], blankVertices.front());
                            std::swap(blankInactive[j], blankInactive.front());
                            ++skipStart;
                        }
                    }
                }
            }
            unsetActive(v);
        }
    }
    return changed;
}

bool PdsState::activateNecessaryNodes() {
    std::vector<Vertex> blankVertices;
    std::vector<Vertex> necessary;
    for (auto v: m_graph.vertices()) {
        if (isBlank(v)) blankVertices.push_back(v);
    }
    for (size_t i = blankVertices.size(); i--;) {
        setActive(blankVertices[i]);
    }
    assert(allObserved());
    for (size_t i = 0; i < blankVertices.size(); ++i) {
        unsetActive(blankVertices[i]);
        if (!allObserved()) {
            necessary.push_back(blankVertices[i]);
        }
        setActive(blankVertices[i]);
    }
    size_t i = 0, j = 0;
    while (i < blankVertices.size() && j < necessary.size()) {
        if (blankVertices[i] == necessary[j]) {
            ++j; ++i;
        } else {
            unsetActive(blankVertices[i]);
            ++i;
        }
    }
    for (; i < blankVertices.size(); ++i) {
        unsetActive(blankVertices[i]);
    }

    return necessary.size() > 0;
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
    auto isObservedEdge = [this](PowerGrid::EdgeDescriptor e) {
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
        m_dependencies.removeEdge(s, t);
        m_dependencies.removeEdge(t, s);
        changed = true;
        auto hasObservedNeighbor = [this](Vertex v) -> bool { return ranges::any_of(m_graph.neighbors(v), [this](auto w) { return isActive(w);});};
        if (!isActive(s)) {
            auto closest = closestActive.at(s);
            if (!hasObservedNeighbor(s)) {
                addEdge(s, closest);
            }
            m_dependencies.addEdge(closest, s);
        }
        if (!isActive(t)) {
            auto closest = closestActive.at(t);
            if (!hasObservedNeighbor(t)) {
                addEdge(t, closest);
            }
            m_dependencies.addEdge(closest, t);
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

SolveState combineSolveState(SolveState first, SolveState second) {
    if (first == SolveState::Infeasible || second == SolveState::Infeasible) {
        return SolveState::Infeasible;
    } else if (first == SolveState::Timeout || second == SolveState::Timeout) {
        return SolveState::Timeout;
    } else if (first == SolveState::Other || second == SolveState::Other) {
        return SolveState::Other;
    } else if (first == SolveState::Heuristic || second == SolveState::Heuristic) {
        return SolveState::Heuristic;
    } else {
        return first;
    }
}

}
