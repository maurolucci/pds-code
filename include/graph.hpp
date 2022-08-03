//
// Created by max on 03.08.22.
//

#ifndef PDS_GRAPH_HPP
#define PDS_GRAPH_HPP

#include <optional>
#include <boost/graph/adjacency_list.hpp>
#include "utility.hpp"

namespace pds {

template<class BoostGraph>
class WrappedGraph {
    BoostGraph m_graph;
    //boost::adjacency_list<OutEdgeList, VertexList, Directed, VertexProperties, EdgeProperties, GraphProperties, EdgeList> m_graph;
public:
    using vertex_descriptor = typename boost::graph_traits<decltype(m_graph)>::vertex_descriptor;
    using edge_descriptor = typename boost::graph_traits<decltype(m_graph)>::edge_descriptor;
    using directed_category = typename boost::graph_traits<decltype(m_graph)>::directed_category;
    using edge_parallel_category = typename boost::graph_traits<decltype(m_graph)>::edge_parallel_category;
    using traversal_category = typename boost::graph_traits<decltype(m_graph)>::traversal_category;

    using out_edge_iterator = typename boost::graph_traits<decltype(m_graph)>::out_edge_iterator;
    using degree_size_type = typename boost::graph_traits<decltype(m_graph)>::degree_size_type;

    WrappedGraph() = default;
    WrappedGraph(const WrappedGraph& other) = default;
    WrappedGraph(WrappedGraph&& other) = default;

    //template<class... T>
    //WrappedGraph(T&&... args) : m_graph(std::forward(args)...) { }

    auto out_edges(vertex_descriptor v) const {
        return pds::range_pair(boost::out_edges(v, m_graph));
    }

    vertex_descriptor source(edge_descriptor edge) const {
        return boost::source(edge, m_graph);
    }

    vertex_descriptor target(edge_descriptor edge) const {
        return boost::target(edge, m_graph);
    }

    static vertex_descriptor null_vertex() {
        return boost::graph_traits<decltype(m_graph)>::null_vertex();
    }

    degree_size_type out_degree(vertex_descriptor v) const {
        return boost::out_degree(v, m_graph);
    }

    using in_edge_iterator = typename decltype(m_graph)::in_edge_iterator;

    auto in_edges(vertex_descriptor v) const {
        return pds::range_pair(boost::in_edges(v, m_graph));
    }

    degree_size_type in_degree(vertex_descriptor v) const {
        return boost::in_degree(v, m_graph);
    }

    degree_size_type degree(vertex_descriptor v) const {
        return boost::degree(v, m_graph);
    }

    using adjacency_iterator = typename decltype(m_graph)::adjacency_iterator;

    auto adjacent_vertices(vertex_descriptor vertex) const {
        return range_pair(boost::adjacent_vertices(vertex, m_graph));
    }

    using vertex_iterator = typename decltype(m_graph)::vertex_iterator;
    using vertices_size_type = typename decltype(m_graph)::vertices_size_type;

    decltype(range_pair(boost::vertices(m_graph))) vertices() const {
        return range_pair(boost::vertices(m_graph));
    }

    vertices_size_type num_vertices() const {
        return boost::num_vertices(m_graph);
    }

    using edge_iterator = typename decltype(m_graph)::edge_iterator;
    using edges_size_type = typename decltype(m_graph)::edges_size_type;

    auto edges() const {
        return range_pair(boost::edges(m_graph));
    }

    edges_size_type num_edges() const {
        return boost::num_edges(m_graph);
    }

    std::optional<edge_descriptor> edge(vertex_descriptor u, vertex_descriptor v) const {
        auto edge = boost::edge(u, v, m_graph);
        if (edge.second) {
            return {edge.first};
        } else {
            return {};
        }
    }

    using vertex_property_type = typename decltype(m_graph)::vertex_property_type;

    template<typename... T>
    vertex_descriptor add_vertex(T&&... args) {
        return boost::add_vertex(vertex_property_type(std::forward(args)...), m_graph);
    }

    void clear_vertex(vertex_descriptor vertex) {
        boost::clear_vertex(vertex, m_graph);
    }

    void remove_vertex(vertex_descriptor vertex) {
        boost::remove_vertex(vertex, m_graph);
    }

    using edge_property_type = typename decltype(m_graph)::edge_property_type;

    template<typename... T>
    std::optional<edge_descriptor> add_edge(vertex_descriptor u, vertex_descriptor v, T&&... args) {
        auto edge = boost::add_edge(u, v, edge_property_type(std::forward(args)...), m_graph);
        if (edge.second) {
            return {edge.first};
        } else {
            return {};
        }
    }

    void remove_edge(vertex_descriptor u, vertex_descriptor v) {
        boost::remove_edge(u, v, m_graph);
    }

    void remove_edge(edge_descriptor e) {
        boost::remove_edge_tag(e, m_graph);
    }

    void remove_edge(out_edge_iterator e) {
        boost::remove_edge(e, m_graph);
    }

    auto& operator[](vertex_descriptor vertex) {
        return m_graph[vertex];
    }

    const auto& operator[](vertex_descriptor vertex) const {
        return m_graph[vertex];
    }

    template<class Property>
    auto getProperty(Property prop) {
        return get(prop, m_graph);
    }

    template<class Property>
    auto getProperty(Property prop) const {
        return get(prop, m_graph);
    }

    template<class Property>
    auto& getVertexProperty(const Property& prop, vertex_descriptor vertex) {
        return get(prop, m_graph, vertex);
    }

    template<class Property>
    const auto& getVertexProperty(const Property& prop, vertex_descriptor vertex) const {
        return get(prop, m_graph, vertex);
    }

    template<class Property>
    auto& getEdgeProperty(const Property& prop, edge_descriptor edge) {
        return get(prop, m_graph, edge);
    }

    template<class Property>
    const auto& getEdgeProperty(const Property& prop, edge_descriptor edge) const {
        return get(prop, m_graph, edge);
    }
};

template<class Prop, class Base>
auto get(Prop prop, const WrappedGraph<Base>& graph) {
    return graph.template getProperty(prop);
}

template<class Prop, class Base>
auto get(Prop prop, WrappedGraph<Base>& graph) {
    return graph.template getProperty(prop);
}

template<class Prop, class Base>
auto get(Prop prop, const WrappedGraph<Base>& graph, typename WrappedGraph<Base>::vertex_descriptor vertex) {
    return graph.template getVertexProperty(prop, vertex);
}

template<class Prop, class Base>
auto get(Prop prop, WrappedGraph<Base>& graph, typename WrappedGraph<Base>::vertex_descriptor vertex) {
    return graph.template getVertexProperty(prop, vertex);
}

template<class Prop, class Base>
auto get(Prop prop, const WrappedGraph<Base>& graph, typename WrappedGraph<Base>::edge_descriptor edge) {
    return graph.template getVertexProperty(prop, edge);
}

template<class Prop, class Base>
auto get(Prop prop, WrappedGraph<Base>& graph, typename WrappedGraph<Base>::edge_descriptor edge) {
    return graph.template getVertexProperty(prop, edge);
}

template<class Prop, class Base, class T>
auto put(Prop prop, WrappedGraph<Base>& graph, typename WrappedGraph<Base>::vertex_descriptor vertex, T&& value) {
    return put(prop, graph.m_graph, vertex, std::forward(value));
}

template<class Prop, class Base, class T>
auto put(Prop prop, WrappedGraph<Base>& graph, typename WrappedGraph<Base>::edge_descriptor edge, T&& value) {
    return put(prop, graph.m_graph, edge, std::forward(value));
}

template<
        class OutEdgeList = boost::listS,
        class VertexList = boost::listS,
        class Directed = boost::directedS,
        class VertexProperties = boost::no_property,
        class EdgeProperties = boost::no_property,
        class GraphProperties = boost::no_property,
        class EdgeList = boost::listS>
using AdjGraph = WrappedGraph<boost::adjacency_list<OutEdgeList, VertexList, Directed, VertexProperties, EdgeProperties, GraphProperties, EdgeList>>;

} //namespace pds

namespace boost {
template<class Base>
inline std::pair<typename pds::WrappedGraph<Base>::out_edge_iterator, typename pds::WrappedGraph<Base>::out_edge_iterator> out_edges(typename pds::WrappedGraph<Base>::vertex_descriptor vertex, const pds::WrappedGraph<Base>& graph) {
    auto rng = graph.out_edges(vertex);
    return std::make_pair(ranges::begin(rng), ranges::end(rng));
}

template<class Base>
inline typename pds::WrappedGraph<Base>::vertex_descriptor source(typename pds::WrappedGraph<Base>::edge_descriptor vertex, const pds::WrappedGraph<Base>& graph) {
    return graph.source(vertex);
}

template<class Base>
inline typename pds::WrappedGraph<Base>::vertex_descriptor target(typename pds::WrappedGraph<Base>::edge_descriptor vertex, const pds::WrappedGraph<Base>& graph) {
    return graph.target(vertex);
}

template<class Base>
inline typename pds::WrappedGraph<Base>::degree_size_type out_degree(typename pds::WrappedGraph<Base>::vertex_descriptor vertex, const pds::WrappedGraph<Base>& graph) {
    return graph.out_degree(vertex);
}

template<class Base>
inline std::pair<typename pds::WrappedGraph<Base>::in_edge_iterator, typename pds::WrappedGraph<Base>::in_edge_iterator> in_edges(typename pds::WrappedGraph<Base>::vertex_descriptor vertex, const pds::WrappedGraph<Base>& graph) {
    auto rng = graph.in_edges(vertex);
    return std::make_pair(ranges::begin(rng), ranges::end(rng));
}

template<class Base>
inline typename pds::WrappedGraph<Base>::degree_size_type degree(typename pds::WrappedGraph<Base>::vertex_descriptor vertex, const pds::WrappedGraph<Base>& graph) {
    return graph.degree(vertex);
}

template<class Base>
inline typename pds::WrappedGraph<Base>::degree_size_type in_degree(typename pds::WrappedGraph<Base>::vertex_descriptor vertex, const pds::WrappedGraph<Base>& graph) {
    return graph.in_degree(vertex);
}

template<class Base>
inline std::pair<typename pds::WrappedGraph<Base>::adjacency_iterator, typename pds::WrappedGraph<Base>::adjacency_iterator> adjacent_vertices(typename pds::WrappedGraph<Base>::vertex_descriptor vertex, const pds::WrappedGraph<Base>& graph) {
    auto rng = graph.adjacent_vertices(vertex);
    return std::make_pair(ranges::begin(rng), ranges::end(rng));
}

template<class Base>
inline std::pair<typename pds::WrappedGraph<Base>::vertex_iterator, typename pds::WrappedGraph<Base>::vertex_iterator> vertices(const pds::WrappedGraph<Base>& graph) {
    auto rng = graph.vertices();
    return std::make_pair(ranges::begin(rng), ranges::end(rng));
}

template<class Base>
inline typename pds::WrappedGraph<Base>::vertices_size_type num_vertices(const pds::WrappedGraph<Base>& graph) {
    return graph.num_vertices();
}

template<class Base>
inline std::pair<typename pds::WrappedGraph<Base>::vertex_iterator, typename pds::WrappedGraph<Base>::vertex_iterator> edges(const pds::WrappedGraph<Base>& graph) {
    auto rng = graph.edges();
    return std::make_pair(ranges::begin(rng), ranges::end(rng));
}

template<class Base>
inline typename pds::WrappedGraph<Base>::edges_size_type num_edges(const pds::WrappedGraph<Base>& graph) {
    return graph.num_edges();
}

template<class Base>
inline std::pair<typename pds::WrappedGraph<Base>::edge_descriptor, bool> edge(typename pds::WrappedGraph<Base>::vertex_descriptor u, typename pds::WrappedGraph<Base>::vertex_descriptor v, const pds::WrappedGraph<Base>& graph) {
    auto edge = graph.edge(u, v);
    return std::make_pair(edge.template value_or(typename pds::WrappedGraph<Base>::edge_descriptor{}), edge.has_value());
}

template<class Base>
inline typename pds::WrappedGraph<Base>::vertex_descriptor add_vertex(pds::WrappedGraph<Base>& graph) {
    return graph.template add_vertex();
}

template<class Base>
inline void clear_vertex(typename pds::WrappedGraph<Base>::vertex_descriptor vertex, const pds::WrappedGraph<Base>& graph) {
    return graph.template clear_vertex(vertex);
}

template<class Base>
inline std::pair<typename pds::WrappedGraph<Base>::edge_descriptor, bool> add_edge(typename pds::WrappedGraph<Base>::vertex_descriptor u, typename pds::WrappedGraph<Base>::vertex_descriptor v, pds::WrappedGraph<Base>& graph) {
    auto edge = graph.template add_edge(u, v);
    return std::make_pair(edge.template value_or(typename pds::WrappedGraph<Base>::edge_descriptor{}), edge.has_value());
}

template<class Base>
inline void remove_edge(typename pds::WrappedGraph<Base>::vertex_descriptor u, typename pds::WrappedGraph<Base>::vertex_descriptor v, const pds::WrappedGraph<Base>& graph) {
    return graph.template remove_edge(u, v);
}

template<class Base>
inline void remove_edge(typename pds::WrappedGraph<Base>::edge_descriptor e, const pds::WrappedGraph<Base>& graph) {
    return graph.template remove_edge(e);
}

template<class Base>
inline void remove_edge(typename pds::WrappedGraph<Base>::out_edge_iterator e, const pds::WrappedGraph<Base>& graph) {
    return graph.template remove_edge(e);
}

template<class Base>
inline typename pds::WrappedGraph<Base>::vertex_descriptor add_vertex(const typename pds::WrappedGraph<Base>::vertex_property_type& data, const pds::WrappedGraph<Base>& graph) {
    return graph.template add_vertex(data);
}

template<class Base>
inline std::pair<typename pds::WrappedGraph<Base>::edge_descriptor, bool> add_edge(typename pds::WrappedGraph<Base>::vertex_descriptor u, typename pds::WrappedGraph<Base>::vertex_descriptor v, const typename pds::WrappedGraph<Base>::edge_property_type& data, const pds::WrappedGraph<Base>& graph) {
    auto edge = graph.template add_edge(u, v, data);
    if (edge) {
        std::make_pair(edge.get(), true);
    } else {
        std::make_pair(typename pds::WrappedGraph<Base>::edge_descriptor{}, false);
    }
}

} // namespace boost

#endif //PDS_GRAPH_HPP
