//
// Created by max on 19.07.22.
//

#ifndef PDS_PDS_HPP
#define PDS_PDS_HPP

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/property_map/dynamic_property_map.hpp>

#include <fstream>

#include "utility.hpp"

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
    std::string id;
    bool zero_injection;
};

using PowerGrid = boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, Bus>;

PowerGrid read_graphml(const std::string& filename, bool all_zero_injection = false) {
    PowerGrid graph;
    boost::dynamic_properties attr(boost::ignore_other_properties);
    std::vector<long> long_zi(num_vertices(graph));
    //typename boost::property_map<PowerGrid, boost::vertex_index_t>::type index = get(boost::vertex_index, graph); //
    boost::associative_property_map<std::map<PowerGrid::vertex_descriptor, long>> zero_injection;
    attr.property("zero_injection", zero_injection);//make_vector_property_map(long_zi, graph));
    std::ifstream graph_in(filename);
    boost::read_graphml(graph_in, graph, attr);
    graph_in.close();
    for (auto it = vertices(graph).first; it != vertices(graph).second; ++it) {
        auto v = *it;
        graph[v].zero_injection = true;//zero_injection[v];
    }
    return graph;
}

template<typename T> struct PrintType;
template<class Graph, class ActivePropertyMap, class ObservedPropertyMap>
bool dominate(const Graph& graph, const ActivePropertyMap& active, ObservedPropertyMap& observed) {
    typename boost::property_map<Graph, boost::vertex_index_t>::type index = get(boost::vertex_index, graph);
    std::for_each(vertices(graph).first, vertices(graph).second, [&] (const auto& v) {
        auto x = index[v];

        if (active[v] == PmuState::Active) {
            observed[v] = true;
            auto adjacent = boost::adjacent_vertices(v, graph);
            for (auto it = adjacent.first; it != adjacent.second; ++it) {
                auto windex = index[*it];
                observed[*it] = true;
            }
        }
    });
    return true;
}

template<class Graph, class ObservedPropertyMap>
std::vector<typename Graph::vertex_descriptor> propagation_targets_k(
        Graph& graph,
        const typename Graph::vertex_descriptor& vertex,
        ObservedPropertyMap& observed,
        size_t max_unobserved = 1)
{
    if (!observed[vertex]) return {};
    std::vector<typename Graph::vertex_descriptor> targets;
    for (const auto& w: range_pair(adjacent_vertices(vertex, graph))) {
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

template<class Graph, class ObservedPropertyMap>
bool propagate(const Graph& graph, ObservedPropertyMap& observed, size_t max_unobserved = 1) {
    typename boost::property_map<Graph, boost::vertex_index_t>::type index = get(boost::vertex_index, graph);
    boost::vector_property_map<size_t, decltype(index)> unobserved_degree(num_vertices(graph), index);
    std::vector<typename Graph::vertex_descriptor> queue;
    for (const auto& v: range_pair(vertices(graph))) {
        if (!observed[v]) {
            for (const auto& w: range_pair(adjacent_vertices(v, graph))) {
                unobserved_degree[w] += 1;
            }
        }
    }
    for (const auto& v: range_pair(vertices(graph))) {
        if (observed[v] && unobserved_degree[v] > 0 && unobserved_degree[v] <= max_unobserved) {
            queue.push_back(v);
        }
    }
    while (!queue.empty()) {
        auto v = queue.back();
        queue.pop_back();
        for (const auto& w: range_pair(adjacent_vertices(v, graph))) {
            if (!observed[w]) {
                observed[w] = true;
                if (observed[w] && unobserved_degree[w] > 0 && unobserved_degree[w] <= max_unobserved) {
                    queue.push_back(w);
                }
                for (const auto& u: range_pair(adjacent_vertices(w, graph))) {
                    unobserved_degree[u] -= 1;
                    auto deg = unobserved_degree[u];
                    if (observed[u] && deg > 0 && deg == max_unobserved) {
                        queue.push_back(u);
                    }
                }
            }
        }
    }
    return true;
}

template<class Graph, class ObservedPropertyMap>
bool observed(const Graph& graph, const ObservedPropertyMap& observed) {
    return std::all_of(boost::vertices(graph).first, boost::vertices(graph).second, [&] (const auto& v) {
        return observed[v];
    });
}
}

#endif //PDS_PDS_HPP
