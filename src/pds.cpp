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
    //typename boost::property_map<PowerGrid, boost::vertex_index_t>::type index = get(boost::vertex_index, graph); //
    std::map<PowerGrid::vertex_descriptor, long> zero_injection_data;
    boost::associative_property_map zero_injection(zero_injection_data);
    //auto index_map = boost::get(boost::vertex_index, graph);
    //std::vector<long> long_zi;
    //auto zero_injection = boost::make_iterator_property_map(long_zi.begin())
    auto id_map = [&graph](const PowerGrid::vertex_descriptor& vertex) -> long& { return graph[vertex].id;};
    auto name_map = [&graph](PowerGrid::vertex_descriptor vertex) -> std::string& { return graph[vertex].name; };
    boost::function_property_map<decltype(id_map), PowerGrid::vertex_descriptor> id(id_map);
    boost::function_property_map<decltype(name_map), PowerGrid::vertex_descriptor> name(name_map);
    //auto id_property = graph.getProperty(&Bus::id);
    //auto name_property = graph.getProperty(&Bus::name);
    attr.property("zero_injection", zero_injection);//make_vector_property_map(long_zi, graph));
    attr.property("name", name);
    attr.property("id", id);
    std::ifstream graph_in(filename);
    pds::read_graphml(graph_in, graph, attr);
    graph_in.close();
    for (auto v: graph.vertices()) {
        graph[v].zero_injection = zero_injection[v] || all_zero_injection;
        std::cout << "node id=" << graph[v].name << std::endl;
    }
    return graph;
}

}