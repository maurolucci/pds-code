//
// Created by max on 25.07.22.
//

#ifndef PDS_GUROBI_SOLVE_HPP
#define PDS_GUROBI_SOLVE_HPP

#include <range/v3/all.hpp>

#include "pds.hpp"

namespace pds {
bool solve_pds(const PowerGrid &graph, map<PowerGrid::vertex_descriptor, PmuState> &active, bool output = false);
} //namespace pds
#endif //PDS_GUROBI_SOLVE_HPP
