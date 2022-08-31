//
// Created by max on 31.08.22.
//

#include <setgraph/graph.hpp>
#include <pds.hpp>
#include <fstream>
#include <filesystem>

#include <range/v3/all.hpp>

namespace pds {
inline PowerGrid readEdgeList(const std::string& filename) {
    namespace fs = std::filesystem;
    if (!fs::exists(fs::path(filename))) {
        throw std::runtime_error("file not found: " + filename);
    }
    std::ifstream infile(filename);
    PowerGrid graph;
    while (infile) {
        PowerGrid::vertex_descriptor i, j;
        infile >> i >> j;
        graph.getOrAddVertex(i);
        graph.getOrAddVertex(j);
        graph.addEdge(i, j);
    }
    return graph;
}

inline PowerGrid readPtxt(const std::string& filename) {
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
                    graph.getOrAddVertex(currentNode);
                    graph.getOrAddVertex(target);
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
}

int main(int argc, const char** argv) {
    if (argc != 2) {
        throw std::runtime_error("expected file name");
    }
    std::string filename = argv[1];
    namespace fs = std::filesystem;
    fmt::print("opening {}\n", fs::absolute(filename).string());
    pds::PowerGrid graph;
    if (filename.ends_with(".graph")) {
        graph = pds::readEdgeList(filename);
    } else if (filename.ends_with(".ptxt")) {
        graph = pds::readPtxt(filename);
    } else {
        throw std::runtime_error("unsupported file: " + filename);
    }
    fmt::print("graph n={}, m={}\n", graph.numVertices(), graph.numEdges());
}
