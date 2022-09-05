//
// Created by max on 31.08.22.
//

#ifndef PDS_GRAPHIO_HPP
#define PDS_GRAPHIO_HPP

#include <setgraph/graph.hpp>
#include <pds.hpp>
#include <fstream>
#include <filesystem>
#include <range/v3/all.hpp>

namespace pds {
PowerGrid readEdgeList(const std::string& filename, bool allZeroInjection = false);

PowerGrid readPtxt(const std::string& filename, bool allZeroInjection = false);

PowerGrid readGraphML(const std::string& filename, bool allZeroInjection = false);

PowerGrid readAutoGraph(const std::string& filename, bool allZeroInjection = false);

}

#endif //PDS_GRAPHIO_HPP
