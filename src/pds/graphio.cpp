#include "graphio.hpp"

#include "setgraph/boost_adapter.hpp"
#include <boost/property_map/function_property_map.hpp>
#include "graphml/graphml.hpp"

namespace pds {

PowerGrid readEdgeList(const std::string &filename, bool allZeroInjection) {
    namespace fs = std::filesystem;
    if (!fs::exists(fs::path(filename))) {
        throw std::runtime_error("file not found: " + filename);
    }
    std::ifstream infile(filename);
    PowerGrid graph;
    while (infile) {
        PowerGrid::vertex_descriptor i, j;
        infile >> i >> j;
        graph.getOrAddVertex(i, Bus{.name=std::to_string(i), .id=long(i), .zero_injection=allZeroInjection});
        graph.getOrAddVertex(j, Bus{.name=std::to_string(j), .id=long(j), .zero_injection=allZeroInjection});
        graph.addEdge(i, j);
    }
    return graph;
}

PowerGrid readPtxt(const std::string &filename, bool allZeroInjection) {
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
                    graph.getOrAddVertex(currentNode,
                                         Bus{.name=std::to_string(currentNode), .id=long(currentNode), .zero_injection=allZeroInjection});
                    graph.getOrAddVertex(target,
                                         Bus{.name=std::to_string(target), .id=long(target), .zero_injection=allZeroInjection});
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

PowerGrid readGraphML(const std::string &filename, bool all_zero_injection) {
    PowerGrid graph;
    boost::dynamic_properties attr(boost::ignore_other_properties);
    map<PowerGrid::vertex_descriptor, long> zero_injection_data;
    boost::associative_property_map zero_injection(zero_injection_data);
    auto id_map = [&graph](const PowerGrid::vertex_descriptor &vertex) -> long & { return graph[vertex].id; };
    auto name_map = [&graph](PowerGrid::vertex_descriptor vertex) -> std::string & { return graph[vertex].name; };
    boost::function_property_map<decltype(id_map), PowerGrid::vertex_descriptor> id(id_map);
    boost::function_property_map<decltype(name_map), PowerGrid::vertex_descriptor> name(name_map);
    attr.property("zero_injection", zero_injection);//make_vector_property_map(long_zi, graph));
    attr.property("name", name);
    attr.property("id", id);
    std::ifstream graph_in(filename);
    read_graphml(graph_in, graph, attr);
    graph_in.close();
    for (auto v: graph.vertices()) {
        graph[v].zero_injection = zero_injection[v] || all_zero_injection;
    }
    return graph;
}

struct ParseError : std::exception {
private:
    std::string reason;
public:
    ParseError(const std::string& reason, size_t line) : reason(fmt::format("{}: {}", line, reason)) { }
    ParseError(const ParseError&) = default;
    ParseError(ParseError&&) = default;


    const char * what() const noexcept override {
        return reason.c_str();
    }
};

std::vector<std::string_view> split(const std::string_view& in, char delim) {
    std::vector<std::string_view> pieces;
    size_t start = 0;
    size_t end;
    while (start != std::string::npos) {
        end = in.find(delim, start);
        pieces.emplace_back(in.substr(start, end - start));
        if (end == std::string::npos) {
            start = std::string::npos;
        } else {
            start = end + 1;
        }

    }
    return pieces;
}

int parseInt(const std::string_view& token, size_t line) {
    try {
        return std::stoi(std::string{token});
    } catch (std::invalid_argument) {
        throw ParseError { fmt::format("not a number: {}", token), line };
    } catch (std::out_of_range) {
        throw ParseError { fmt::format("invalid number: {}", token), line };
    }
}

PowerGrid readPds(const std::string& filename, bool allZeroInjection) {
    std::ifstream infile(filename);
    std::string line;
    size_t lineno = 0;
    std::getline(infile, line);
    ++lineno;
    if (!line.starts_with("PDS Instance V1.0")) {
        throw ParseError{"Invalid Header", 1};
    }
    std::getline(infile, line);
    ++lineno;
    size_t numVertices, numEdges;
    if (!line.starts_with("G")) {
        throw ParseError{"Expected graph information", 2};
    } else {
        auto pieces = split(line, ' ');
        if (pieces.size() != 3) throw ParseError{"incomplete graph information", lineno};
        if (pieces[0] != "G") throw ParseError{"unexpected line", lineno};
        numVertices = parseInt(pieces[1], lineno);
        numEdges = parseInt(pieces[2], lineno);
    }

    PowerGrid graph;
    while (infile) {
        std::getline(infile, line);
        ++lineno;
        auto pieces = split(line, ' ');
        if (pieces.size() > 0) {
            if (pieces[0] == "V") {
                if (pieces.size() < 3) throw ParseError{"too little data", lineno};
                int v = parseInt(pieces[1], lineno);
                if (graph.hasVertex(v)) throw ParseError{"duplicate vertex definition", lineno};
                auto &data = graph.getOrAddVertex(v);
                data.name = pieces[1];
                data.zero_injection = allZeroInjection;
                data.pmu = PmuState::Blank;
                if (pieces.size() == 4) {
                    for (auto c: pieces[3]) {
                        switch (c) {
                            case 'Z':
                                data.zero_injection = true;
                                break;
                            case 'A':
                                data.pmu = PmuState::Active;
                                break;
                            case 'I':
                                data.pmu = PmuState::Inactive;
                                break;
                            default:
                                throw ParseError{"unknown vertex type", lineno};
                        }
                    }
                } else if (pieces.size() > 4) {
                    throw ParseError{"too many vertex arguments", lineno};
                }
            } else if (pieces[0] == "E") {
                if (pieces.size() != 3) {
                    throw ParseError{"", lineno};
                }
                auto s = parseInt(pieces[1], lineno);
                auto t = parseInt(pieces[2], lineno);
                if (!graph.hasVertex(s)) throw ParseError { fmt::format("undefined start: {}", s), lineno};
                if (!graph.hasVertex(t)) throw ParseError { fmt::format("undefined end: {}", s), lineno};
                graph.addEdge(s, t);
            } else if(pieces[0] != "") {
                throw ParseError{fmt::format("unexpected line: {}", line), lineno};
            }
        }
    }
    if (graph.numVertices() != numVertices) throw ParseError { fmt::format("expected {} vertices but got {}", numVertices, graph.numVertices()), lineno};
    if (graph.numEdges() != numEdges) throw ParseError { fmt::format("expected {} vertices but got {}", numVertices, graph.numVertices()), lineno};
    return graph;
}

void writePds(const PowerGrid& grid, const std::string& filename) {
    FILE* outfile = fopen(filename.c_str(), "wb");
    fmt::print(outfile, "PDS Instance V1.0\n");
    fmt::print(outfile, "G {} {}\n", grid.numVertices(), grid.numEdges());
    for (auto v: grid.vertices()) {
        const auto& data = grid[v];
        fmt::print(outfile, "V {} \"{}\" ", v, data.name);
        if (data.zero_injection) {
            fmt::print(outfile, "Z");
        }
        switch (data.pmu) {
            case PmuState::Active:
                fmt::print(outfile, "A");
                break;
            case PmuState::Inactive:
                fmt::print(outfile, "I");
                break;
            default:
                break;
        }
        fmt::print(outfile, "\n");
    }
    for (auto [s, t]: grid.edges()) {
        fmt::print(outfile, "E {} {}\n", s, t);
    }
    fclose(outfile);
}

PowerGrid readAutoGraph(const std::string &filename, bool allZeroInjection) {
    if (filename.ends_with(".graphml")) {
        return readGraphML(filename, allZeroInjection);
    } else if (filename.ends_with(".graph")) {
        return readEdgeList(filename, allZeroInjection);
    } else if (filename.ends_with(".ptxt")) {
        return readPtxt(filename, allZeroInjection);
    } else if (filename.ends_with(".pds")) {
        return readPds(filename, allZeroInjection);
    } else {
        throw std::runtime_error("unsupported format: " + filename);
    }
}

} //namespace pds