#ifndef PDS_GUROBI_COMMON_HPP
#define PDS_GUROBI_COMMON_HPP

#include "pds.hpp"

struct GRBModel;
struct GRBVar;
struct GRBEnv;

namespace pds {

struct MIPModel {
    std::unique_ptr<GRBModel> model;
    map<PowerGrid::VertexDescriptor, GRBVar> xi;
    MIPModel();
    MIPModel(MIPModel&& other) = default;
    virtual ~MIPModel();
};

void preloadMIPSolver();
GRBEnv& getEnv();

void relaxMIPModel(MIPModel&);

}

#endif //PDS_GUROBI_COMMON_HPP
