#include "graphio.hpp"

#include "setgraph/boost_adapter.hpp"
#include <boost/property_map/function_property_map.hpp>
#include "graphml/graphml.hpp"

pds::PowerGrid pds::readEdgeList(const std::string &filename, bool allZeroInjection) {
    namespace fs = std::filesystem;
    if (!fs::exists(fs::path(filename))) {
        throw std::runtime_error("file not found: " + filename);
    }
    std::ifstream infile(filename);
    PowerGrid graph;
    while (infile) {
        PowerGrid::vertex_descriptor i, j;
        infile >> i >> j;
        graph.getOrAddVertex(i, Bus {.name=std::to_string(i), .id=long(i), .zero_injection=allZeroInjection});
        graph.getOrAddVertex(j, Bus {.name=std::to_string(j), .id=long(j), .zero_injection=allZeroInjection});
        graph.addEdge(i, j);
    }
    return graph;
}

pds::PowerGrid pds::readPtxt(const std::string &filename, bool allZeroInjection) {
    namespace fs = std::filesystem;
    if (!fs::exists(fs::path(filename))) {
        throw std::runtime_error("file not found: " + filename);
    }
    std::ifstream infile(filename);
    PowerGrid graph;
    enum class State {
        Init,
        Graph,
        Node,
        Edges,
        ExpectList,
        IncidenceList,
    };
    State state = State::Init;
    std::string line;
    size_t currentNode = -1;
    size_t expectedEdges = 0;
    using namespace std::string_literals;
    while (std::getline(infile, line)) {
        std::string lower = line;
        ranges::transform(lower, lower.begin(), ::tolower);
        if (lower.starts_with("number edges containting node")) { //sic
            state = State::Node;
            expectedEdges = 0;
            std::string nodeStr = lower.substr("number edges containting node"s.size());
            currentNode = std::stoi(nodeStr);
        } else if (lower.starts_with("number of vertices")) {
            if (state != State::Init) throw std::runtime_error("unexpected line: " + line);
            state = State::Graph;
        } else if (lower.starts_with("edges")) {
            if (state != State::Graph) throw std::runtime_error("unexpected line: " + line);
            state = State::Graph;
        } else if (lower.starts_with("list")) {
            if (state != State::ExpectList) throw std::runtime_error("unexpected line: " + line);
            state = State::IncidenceList;
        } else {
            if (state == State::Node) {
                trim(lower);
                expectedEdges = std::stoi(lower);
                state = State::ExpectList;
            } else if (state == State::IncidenceList) {
                if (expectedEdges) {
                    trim(lower);
                    size_t target = std::stoi(lower);
                    graph.getOrAddVertex(currentNode, Bus {.name=std::to_string(currentNode), .id=long(currentNode), .zero_injection=allZeroInjection});
                    graph.getOrAddVertex(target, Bus {.name=std::to_string(target), .id=long(target), .zero_injection=allZeroInjection});
                    graph.addEdge(currentNode, target);
                    --expectedEdges;
                } else {
                    throw std::runtime_error("unexpected line: " + line);
                }
            } else if (state == State::Graph) {
                //trim(lower);
                //size_t numVertices = std::stoi(lower);
                //unused(numVertices);
                // ignore
            } else {
                throw std::runtime_error("unexpected line: " + line);
            }
        }
    }
    return graph;
}

pds::PowerGrid pds::readGraphML(const std::string &filename, bool all_zero_injection) {
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

pds::PowerGrid pds::readAutoGraph(const std::string &filename, bool allZeroInjection) {
    if (filename.ends_with(".graphml")) {
        return readGraphML(filename, allZeroInjection);
    } else if (filename.ends_with(".graph")) {
        return readEdgeList(filename, allZeroInjection);
    } else if (filename.ends_with(".ptxt")) {
        return readPtxt(filename, allZeroInjection);
    } else {
        throw std::runtime_error("unsupported format: " + filename);
    }
}
