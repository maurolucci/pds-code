#include "graphio.hpp"

#include <tinyxml2.h>

namespace pds {

PowerGrid readEdgeList(const std::string &filename, bool allZeroInjection) {
    namespace fs = std::filesystem;
    if (!fs::exists(fs::path(filename))) {
        throw std::runtime_error("file not found: " + filename);
    }
    std::ifstream infile(filename);
    PowerGrid graph;
    map<size_t, PowerGrid::VertexDescriptor> vertices;
    auto getVertex = [&graph, &vertices,allZeroInjection](size_t v) {
        if (vertices.contains(v)) { return vertices[v]; }
        vertices[v] = graph.addVertex(Bus{.name=std::to_string(v), .id=static_cast<long>(v), .zero_injection=allZeroInjection});
        return vertices[v];
    };
    while (infile) {
        PowerGrid::VertexDescriptor i, j;
        infile >> i >> j;
        graph.addEdge(getVertex(i), getVertex(j));
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
    map<size_t, PowerGrid::VertexDescriptor> vertices;
    auto getVertex = [&graph, &vertices,allZeroInjection](size_t v) {
        if (vertices.contains(v)) { return vertices[v]; }
        vertices[v] = graph.addVertex(Bus{.name=std::to_string(v), .id=static_cast<long>(v), .zero_injection=allZeroInjection});
        return vertices[v];
    };
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
                    graph.addEdge(getVertex(currentNode), getVertex(target));
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

namespace {
struct GraphMLAttribute {
    std::string type;
    std::string name;
};
bool parseBool(const std::string& text) {
    auto boolString = ranges::transform_view(text, [](auto c) { return std::tolower((unsigned char) c); }) | ranges::to<std::string>;
    if (boolString == "false" || boolString == "0" || boolString == "" || boolString == "no") return false;
    return true;
}
}

PowerGrid readGraphML(const std::string &filename, bool all_zero_injection) {
    using namespace std::string_literals;
    using std::string;
    PowerGrid outgraph;
    tinyxml2::XMLDocument doc;
    auto err = doc.LoadFile(filename.c_str());
    if (err != tinyxml2::XMLError::XML_SUCCESS) {
        throw ParseError(doc.ErrorIDToName(err), 0);
    }
    auto graphml = doc.FirstChildElement("graphml");
    if (!graphml) throw ParseError("no graphml content", doc.GetLineNum());
    auto graph = graphml->FirstChildElement("graph");
    if (!graph) throw ParseError("no graph", graphml->GetLineNum());
    auto attributes = graphml->FirstChildElement("key");
    map<string, GraphMLAttribute> nodeAttributes;
    while (attributes) {
        const char *id, *type, *element, *name;
        auto queryAttribute = [&doc, &attributes](const char* name, const char** attr) -> bool {
            auto err = attributes->QueryAttribute(name, attr);
            if (err) {
                fmt::print(stderr, "[WARN] cannot read attribute {}: {}", name, doc.ErrorIDToName(err));
                return false;
            } else { return true; }
        };
        if (queryAttribute("id", &id)
            && queryAttribute("for", &element)
            && queryAttribute("attr.type", &type)
            && queryAttribute("attr.name", &name)
        ) {
            nodeAttributes[id] = {type, name};
        }
        attributes = attributes->NextSiblingElement("key");

    }
    auto node = graph->FirstChildElement("node");
    map<std::string, PowerGrid::VertexDescriptor> vertices;
    while (node) {
        const char* key;
        if (node->QueryAttribute("id", &key)) throw ParseError("invalid node", node->GetLineNum());
        auto vertex = outgraph.addVertex(Bus {
            .name=key,
            .id=static_cast<long>(vertices.size()),
            .zero_injection=all_zero_injection,
            .pmu=PmuState::Blank
        });
        vertices[key] = vertex;
        auto data = node->FirstChildElement("data");
        while (data) {
            auto text = data->GetText();
            if (!text) throw ParseError("could not read text", data->GetLineNum());
            if (!data->QueryAttribute("key", &key)) {
                if ("zero_injection"s == nodeAttributes[key].name) {
                    outgraph.getVertex(vertex).zero_injection = all_zero_injection || parseBool(text);
                } else if ("name"s == nodeAttributes[key].name) {
                    outgraph.getVertex(vertex).name = text;
                }
            }
            data = data->NextSiblingElement("data");
        }
        node = node->NextSiblingElement("node");
    }
    auto edge = graph->FirstChildElement("edge");
    while (edge) {
        const char* source, * target;
        if (edge->QueryAttribute("source", &source)) throw ParseError("cannot parse edge source", edge->GetLineNum());
        if (edge->QueryAttribute("target", &target)) throw ParseError("cannot parse edge target", edge->GetLineNum());
        if (!vertices.contains(source) || !vertices.contains(target)) throw ParseError("invalid edge", edge->GetLineNum());
        outgraph.addEdge(vertices[source], vertices[target]);
        edge = edge->NextSiblingElement("edge");
    }
    return outgraph;
}

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
    map<int, PowerGrid::VertexDescriptor> vertices;
    while (infile) {
        std::getline(infile, line);
        ++lineno;
        auto pieces = split(line, ' ');
        if (pieces.size() > 0) {
            if (pieces[0] == "V") {
                if (pieces.size() < 3) throw ParseError{"too little data", lineno};
                int v = parseInt(pieces[1], lineno);
                if (vertices.contains(v)) throw ParseError{"duplicate vertex definition", lineno};
                vertices[v] = graph.addVertex();
                auto &data = graph.getVertex(vertices[v]);
                data.id = v;
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
                if (!vertices.contains(s)) throw ParseError { fmt::format("undefined start: {}", s), lineno};
                if (!vertices.contains(t)) throw ParseError { fmt::format("undefined end: {}", s), lineno};
                graph.addEdge(vertices[s],  vertices[t]);
            } else if(pieces[0] != "") {
                throw ParseError{fmt::format("unexpected line: {}", line), lineno};
            }
        }
    }
    if (graph.numVertices() != numVertices) throw ParseError { fmt::format("expected {} vertices but got {}", numVertices, graph.numVertices()), lineno};
    if (graph.numEdges() != numEdges) throw ParseError { fmt::format("expected {} edges but got {}", numEdges, graph.numEdges()), lineno};
    return graph;
}

void writePds(const PowerGrid& grid, const std::string& filename) {
    FILE* outfile = fopen(filename.c_str(), "wb");
    fmt::print(outfile, "PDS Instance V1.0\n");
    fmt::print(outfile, "G {} {}\n", grid.numVertices(), grid.numEdges());
    assert(grid.numVertices() == size_t(ranges::distance(grid.vertices())));
    for (auto v: grid.vertices()) {
        const auto& data = grid[v];
        fmt::print(outfile, "V {} \"{}\" ", data.id, data.name);
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
    auto id = [&grid](auto v) {return grid[v].id; };
    for (auto [s, t]: grid.edges() | ranges::views::transform([&grid](auto e) { return grid.endpoints(e); })) {
        fmt::print(outfile, "E {} {}\n", id(s), id(t));
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