//
// Created by max on 31.08.22.
//

#include "graphio.hpp"

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
