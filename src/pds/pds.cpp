//
// Created by max on 01.08.22.
//
#include "pds.hpp"

#include "utility.hpp"

namespace pds {

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
                // cannot call setActive on already active vertex
                // setBlank on active vertex requires previous call to setActive
                m_graph[v].pmu = PmuState::Blank;
                setActive(v);
                break;
            case PmuState::Inactive:
                m_graph[v].pmu = PmuState::Blank;
                setInactive(v);
                break;
            default:
                break;
        }
    }
    assert(ranges::distance(m_graph.vertices() | ranges::views::filter([this](auto v) { return isActive(v); })) == (ssize_t)m_numActive);
    assert(ranges::distance(m_graph.vertices() | ranges::views::filter([this](auto v) { return isInactive(v); })) == (ssize_t)m_numInactive);
    assert(ranges::all_of(m_graph.vertices(), [this](auto v) { return m_unobserved_degree.at(v) == ranges::distance(m_graph.neighbors(v) | ranges::views::filter([this](auto v) { return !isObserved(v);})); }));
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
        assert(isActive(vertex) || !isObserved(vertex) || m_dependencies.inDegree(vertex) > 0);
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

bool PdsState::disableLowDegreeRecursive(PdsState::Vertex start, VertexSet &seen) {
    bool changed = false;
    seen.insert(start);
    for (auto w: m_graph.neighbors(start)) {
        if (!seen.contains(w)) {
            changed |= disableLowDegreeRecursive(w, seen);
        }
    }
    auto hasBlankNeighbor = [this] (auto v) {
        return ranges::any_of(m_graph.neighbors(v), [this](auto w) { return isActive(w) || isBlank(w);} );
    };
    if (m_graph.degree(start) <= 2 && hasBlankNeighbor(start) && isZeroInjection(start) && isBlank(start)) {
        setInactive(start);
        changed = true;
    }
    return changed;
}

bool PdsState::setBlank(PdsState::Vertex vertex) {
    if (isInactive(vertex)) {
        --m_numInactive;
        m_graph[vertex].pmu = PmuState::Blank;
    } else if (isActive(vertex)) {
        unsetActive(vertex);
    }
    return allObserved();
}

bool PdsState::setActive(PdsState::Vertex vertex) {
    if (!isActive(vertex)) {
        if (isInactive(vertex)) {
            --m_numInactive;
        }
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
                    if (m_dependencies.hasEdge(v, w)) {
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
    if (!isInactive(vertex)) {
        if (isActive(vertex)) unsetActive(vertex);
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
            Vertex neighbor{};
            for (auto w: m_graph.neighbors(v)) {
                neighbor = w; break;
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
    //set<PdsState::Vertex> seen;
    assert(m_seen.empty());
    for (auto v: m_graph.vertices()) {
        if (m_graph.degree(v) >= 3) {
            changed |= disableLowDegreeRecursive(v, m_seen);
        }
    }
    for (auto v: m_graph.vertices()) {
        m_seen.erase(v);
    }
    assert(m_seen.empty());
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
            if (isBlank(v)) {
                if (isBlank(x) || isBlank(y)) {
                    setInactive(v);
                }
            }
            if (isInactive(v)) {
                if (!m_graph.edge(x, y)) {
                    auto otherNeighbor = [this] (auto i, auto j) {
                        return (m_graph.neighbors(i) | ranges::views::filter([this,j](auto v) { return v != j && !isObserved(v);}) | ranges::to_vector)[0];
                    };
                    if (!isObserved(v) && isObserved(x) && isInactive(x) && isZeroInjection(x) && m_unobserved_degree[x] == 2 && m_graph.degree(otherNeighbor(x, v)) == 2 && isZeroInjection(otherNeighbor(x, v))) {
                        // v must be unobserved if x is observed
                        auto z = otherNeighbor(x, v);
                        if (!m_graph.hasEdge(y, z)) {
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
                        if (!m_graph.hasEdge(x, z)) {
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
                } else {
                    for (auto [s, t]: {std::tie(x, y), std::tie(y, x)}) {
                        if (isBlank(t) && m_graph.degree(s) == 2 && isInactive(s)) {
                            setActive(t);
                            changed = true;
                        }
                    }
                }
                if (!isZeroInjection(x) && !isZeroInjection(y)) {
                    for (auto [s, t]: {std::tie(x, y), std::tie(y, x)}) {
                        if (isInactive(s) && isBlank(t)) {
                            setActive(t);
                            changed = true;
                        }
                    }
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

PdsState::Vertex findActiveParent(const PdsState& state, PdsState::Vertex vertex) {
    while (!state.isActive(vertex)) {
        assert(!ranges::empty(state.observationGraph().inNeighbors(vertex)));
        vertex = *(state.observationGraph().inNeighbors(vertex).begin());
    }
    return vertex;
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
    auto edges = m_graph.edges()
            | ranges::views::filter(isObservedEdge)
            | ranges::views::transform([this](auto e) { return graph().endpoints(e); })
            | ranges::to<std::vector>();
    if (edges.empty()) return false;
    //auto closestActive = findClosestActive(*this);

#ifndef NDEBUG
    for (auto v: graph().vertices()) {
        assert(!isObserved(v) || isActive(v) || !ranges::empty(m_dependencies.inNeighbors(v)));
    }
#endif
    for (auto e: edges) {
        auto [s, t] = e;
        if ((!isActive(s) && isActive(t) && isObservingEdge(t, s)) || (!isActive(t) && isActive(s) && isObservingEdge(s, t))) { continue; }
        m_graph.removeEdge(s, t);
        auto closestS = findActiveParent(*this, s);
        auto closestT = findActiveParent(*this, t);
        m_dependencies.removeEdge(s, t);
        m_dependencies.removeEdge(t, s);
        changed = true;
        auto hasObservedNeighbor = [this](Vertex v) -> bool { return ranges::any_of(m_graph.neighbors(v), [this](auto w) { return isActive(w);});};
        if (!isActive(s)) {
            if (!hasObservedNeighbor(s)) {
                addEdge(s, closestS);
            }
            m_dependencies.addEdge(closestS, s);
        }
        if (!isActive(t)) {
            if (!hasObservedNeighbor(t)) {
                addEdge(t, closestT);
            }
            m_dependencies.addEdge(closestT, t);
        }
        assert(isActive(s) || !ranges::empty(m_dependencies.inNeighbors(s)));
        assert(isActive(t) || !ranges::empty(m_dependencies.inNeighbors(t)));
    }
#ifndef NDEBUG
    for (auto v: graph().vertices()) {
        assert(!isObserved(v) || isActive(v) || !ranges::empty(m_dependencies.inNeighbors(v)));
    }
#endif

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
        set<PdsState::Vertex>& currentComponent
) {
    std::vector<PdsState::Vertex> queue;
    queue.push_back(start);
    seen.insert(start);
    while (!queue.empty()) {
        auto current = queue.back();
        queue.pop_back();
        currentComponent.insert(current);
        for (auto w: state.graph().neighbors(current)) {
            if (!seen.contains(w)) {
                seen.insert(w);
                queue.push_back(w);
            }
        }
    }
}

std::vector<PdsState> PdsState::subproblems() const {
    std::vector<PdsState> components;
    set<Vertex> seen;
    for (auto v: m_graph.vertices()) {
        if (!isActive(v) && !seen.contains(v)) {
            set<Vertex> currentComponent;
            recurseComponent(v, *this, seen, currentComponent);
            PowerGrid subgraph{graph()};
            for (auto x: subgraph.vertices()) {
                if (!currentComponent.contains(x)) {
                    subgraph.removeVertex(x);
                }
            }
            size_t numVertices = subgraph.numVertices(); unused(numVertices);
#ifndef USE_HASHMAP
            subgraph.compress();
#endif
            assert(numVertices == subgraph.numVertices());
            assert(ssize_t(numVertices) == ranges::distance(subgraph.vertices()));
            assert(ranges::all_of(subgraph.vertices(), [&subgraph](auto v) { return subgraph.hasVertex(v);}));
            components.emplace_back(std::move(subgraph));
        }
    }
    return components;
}

void PdsState::applySubsolution(const pds::PdsState &other) {
#ifdef USE_HASHMAP
    for (auto v: other.graph().vertices()) {
        if (other.isActive(v)) {
            setActive(v);
        }
    }
#else
    if (graph().numVertices() == 0 || other.graph().numVertices() == 0) return;
    auto ownVertices = graph().vertices();
    auto ownIt = ranges::begin(ownVertices);
    auto ownEnd = ranges::end(ownVertices);
    auto otherVertices = other.graph().vertices();
    auto otherIt = ranges::begin(otherVertices);
    auto otherEnd = ranges::end(otherVertices);

    while (ownIt != ownEnd && otherIt != otherEnd) {
        const auto ownId = graph().getVertex(*ownIt).id;
        const auto otherId = other.graph().getVertex(*otherIt).id;
        if (ownId < otherId) ++ownIt;
        else if (otherId < ownId) ++otherIt;
        else {
            if (other.isActive(*otherIt)) {
                if (!isInactive(*ownIt)) {
                    setActive(*ownIt);
                } else {
                    fmt::print(stderr, "!!!inactive vertex set to active in subproblem!!!\n");
                }
            }
            ++ownIt; ++otherIt;
        }
    }

#endif
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
