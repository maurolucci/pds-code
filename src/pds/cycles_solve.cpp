#include "cycles_solve.hpp"

#include <gurobi_c++.h>
#include <range/v3/all.hpp>
#include <utility>

#include "pdssolve.hpp"
#include "gurobi_common.hpp"

namespace pds {

// std::vector<VertexList> initializeCycles(PdsState& state, int variant, VertexSet& seen) {
//     switch (variant) {
//         default:
//             return {};
//     }
// }

// auto violatedCycles(PdsState& lastSolution, int variant, double remainingTimeout, VertexSet& seen, int output) -> std::vector<VertexList> {
//     switch(variant) {
//         default:
//             return {};
//     }
// }

struct Callback : public GRBCallback {

    // Lower bound
    // TODO: que es
    size_t* lower;

    // An upper bound for the power dominating number
    // Number of active vertices + undecided vertices
    size_t* upper;

    // Map: Vertex -> GRBVar
    VertexMap<GRBVar>* si;
    // TODO: Resto de las variables de arcos dirigidos?

    // Base integer solution (before re-optimization) 
    // It is not a power dominating set
    const PdsState* base;
    
    // Last integer solution (after re-optimizing) 
    // It may or may not be a power dominating set
    PdsState* solution;

    // Last power dominating set found by gurobi
    PdsState* upperBound;

    // Span of undecided vertices
    std::span<PdsState::Vertex> blank;

    // TODO: ??
    int earlyStop;
    
    // Whether to add intermediate cycles 
    int intermediateCycles;

    // Vector of cycles (list of vertices) 
    std::vector<VertexList>* cycles;
    
    // ???
    VertexSet* seen;

    Callback(size_t& lower, size_t& upper, VertexMap<GRBVar>& si,
             const PdsState& base, PdsState& solution, PdsState & upperBound,
             std::span<PdsState::Vertex> blank, int earlyStop,
             int intermediateCycles, std::vector<VertexList>& cycles, VertexSet& seen)
            : lower(&lower), upper(&upper), si(&si), base(&base), solution(&solution),
              upperBound(&upperBound), blank(blank), earlyStop(earlyStop),
              intermediateCycles(intermediateCycles), cycles(&cycles), seen(&seen)
    { }

    void callback() override {

        // General MIP callback
        // Inside the NoRel Heuristic (when solving the root node and inside every node)
        if (where == GRB_CB_MIP) {
            
            if (getIntInfo(GRB_CB_MIP_SOLCNT) > 0) {
                // At least one integer solution was found
    
                // MIP upper bound (+0.5 to avoid numerical issues)
                auto objVal = static_cast<size_t>(getDoubleInfo(GRB_CB_MIP_OBJBST) + 0.5);
                // MIP lower bound
                auto objBound = static_cast<size_t>(getDoubleInfo(GRB_CB_MIP_OBJBND));

                // If earlyStop = 2, this stops the reoptimization as soon as an invalid solution is found,
                // i.e. a non-power dominaintg set, whose objective value is lower or equal to the upper bound
                // (without finishing proving optimality).
                // TODO: no entiendo: (objBound > *lower || earlyStop > 2)
                if (objVal <= *upper && !solution->allObserved() && earlyStop > 1 && (objBound > *lower || earlyStop > 2)) {
                    abort();
                }
            }
        }

        // MIP solution callback
        // An integer solution was found (it does not necessarily improve the incumbent)
        if (where == GRB_CB_MIPSOL) {

            // MIP upper bound (+0.5 to avoid numerical issues)
            // TODO: está ok sumarle 0.5? No será mucho?
            auto objVal = static_cast<size_t>(getDoubleInfo(GRB_CB_MIPSOL_OBJBST) + 0.5);
            // MIP lower bound
            auto objBound = static_cast<size_t>(getDoubleInfo(GRB_CB_MIPSOL_OBJBND));

            if (objVal <= *upper) {
                // The objective value of the solution is as good as the upper bound for the power dominating number

                // Process the solution
                for (auto v: blank) {

                    if (getSolution(si->at(v)) > 0.5) {
                        // If v is selected, then mark v as active 
                        solution->setActive(v);
                    } else if (base->isBlank(v)) {
                        // if v is not selected and v is undecided in the base solution,
                        // then mark v as undecided
                        solution->setBlank(v);
                    } else {
                        // If is not selected and v is active in the base solution,
                        // then mark v as inactive
                        solution->setInactive(v);
                    }
                }

                fmt::print("feasible solution {} <= {}; {} <= {}; {}\n", *lower, objBound, objVal, *upper, solution->allObserved());
                
                if (!solution->allObserved()) {
                    // The solution is not a power dominating set

                    // // TODO: parametros de violatedCycles ???
                    // if (intermediateCycles >= 0) {
                    //     // Cycles must be added for intermediate solutions
                    //     for (auto& f: violatedCycles(*solution, intermediateCycles, 10, *seen, false)) {
                    //         cycles->push_back(std::move(f));
                    //     }
                    // }

                    // If earlyStop = 2, this stops the reoptimization as soon as an invalid solution is found,
                    // i.e. a non-power dominaintg set, whose objective value is lower or equal to the upper bound
                    // (without finishing proving optimality).
                    // TODO: por que se pide objBound > *lower
                    if (earlyStop > 1 && objBound > *lower) { abort(); }
    
                } else if (objVal < *upper) {
                    // The solution is a power dominating set, 
                    // whose objetive value is better than the upper bound for the power dominating number
                    // i.e. Gurobi improved the upper bound for the power dominating number

                    fmt::print("gurobi H {}\n", objVal);
                    *upper = objVal;
                    *upperBound = *solution;
                }
            }
        }
    }
};

auto now() {
    return std::chrono::high_resolution_clock::now();
}

SolveResult solveCycles(
        PdsState &state,                    // Solution of the last reoptimization
        int output,                         // Wheter to log 
        double timeLimit,                   // Time limit
        int variant,                        // TODO: ??
        int cycleInit,                     // Wheter to initialize cycles 
        int greedyUpper,                    // Whether to apply a greedy heuristic to find a pds
        int earlyStop,                      // Whether to stop early
        callback::CycleCallback cycleCB,  // Callback 
        BoundCallback boundCallback,        // TODO: ??
        int intermediateCycles              // Whether to add cycles for intermediate solutions
) {

    // Initialize lastSolution to be empty
    auto lastSolution = state;

    // Initialize feasibleSolution to be empty
    // Best power dominating set found so far
    // It can be updated by the MIP solution callback or the greedy heuristic
    PdsState feasibleSolution = state;
    
    //writePds(lastSolution.graph(), fmt::format("out/0_input.pds"));

    // Undecided vertices in the solution of the last reoptimization
    auto blankVertices = state.graph().vertices()
                         | ranges::views::filter([&state] (auto v) { return state.isBlank(v); })
                         | ranges::to<std::vector<PdsState::Vertex>>;

    // Set the lower bound to zero
    size_t lowerBound = 0;

    // The set of active or undecided vertices is a trivial power dominating set
    // This set gives an upper bound for the power dominating number
    // The upper bound can be updated by the MIP solution callback, the greedy heuristic, or this function
    size_t upperBound = state.numActive() + state.numBlank();

    if (state.allObserved()) {
        // The solution of the last reoptimization is a power dominating set
        return {state.numActive(), state.numActive(), SolveState::Optimal};
    } else if (state.numBlank() == 0) {
        // The solution of the last reoptimization is not a power dominating set
        // and there is no undecided vertices
        return {state.graph().numVertices(), 0, SolveState::Infeasible};
    }

    try {

        // TODO: ??
        VertexSet seen;

        // Initialize cycles
        // auto cycles = initializeCycles(state, cycleInit, seen);
        std::vector<VertexList> cycles();
        cycleInit = 0;
        variant = 0;

        // Intitialize model
        GRBModel model(getEnv());

        // Set parameters
        model.set(GRB_IntParam_LogToConsole, false);
        model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        model.set(GRB_IntParam_NumericFocus, 1);
        model.set(GRB_DoubleParam_MIPGap, 1e-6);
        model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);

        if (output) {
            // Set file for logging
            model.set(GRB_StringParam_LogFile, "gurobi.log");
        }

        // Add variables
        VertexMap<GRBVar> sv;
        EdgeMap<GRBVar> ye;
        for (auto v: state.graph().vertices()) {

            // Variable: s_v
            sv.emplace(v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("s_{}", v)));
            if (state.isActive(v)) {
                model.addConstr(sv.at(v) == 1);
            }
            if (state.isInactive(v)) {
                model.addConstr(sv.at(v) == 0);
            }

            // Variable: y_e
            if (state.isZeroInjection(v)) {
                for (auto e: state.graph().outEdges(v)) {
                    auto u = state.graph().target(e);
                    ye[v].emplace(u, model.addVar(0.0, 1.0, 0.0, GRB_BINARY, fmt::format("y_{}{}", v, u)));
                }
            }

            // Variable: z_e
            // TODO
        }

        // Add callback
        // TODO: add variables ye and ze
        Callback cb(
                lowerBound, upperBound,
                sv, state, lastSolution, feasibleSolution,
                blankVertices, earlyStop, intermediateCycles, cycles, seen);
        model.setCallback(&cb);

        // Start the clock
        auto startingTime = now();

        // Set the initial status of the model to zero (unsolved)
        int status = 0;

        // Initialize statics for cycles
        // size_t processedCycles = 0;
        size_t totalCycleSize = 0;

        while (true) {
            // Add violated lazy contraints before reoptimization

            // Current time
            auto currentTime = now();

            // Remaining time for timeout
            double remainingTimeout = std::max(0.0, timeLimit - std::chrono::duration_cast<std::chrono::duration<double>>(currentTime - startingTime).count());
            
            // Find cycles
            // TODO :completar
            // auto moreCycles = violatedCycles(lastSolution, variant, remainingTimeout, seen, output);
            std::vector<VertexList> moreCycles();

            // Add the cycles to the set of cycles
            // Recall that the callback may have encountered some cycles in intermediate solutions
            for (auto f: moreCycles) {
                cycles.emplace_back(std::move(f));
            }

            // TODO: creo que nunca se entra a este caso
            // if (forts.empty()) {
            //     // If there is no 
            //     return {state.graph().numVertices(), 0, SolveState::Infeasible};
            // }

            // if (forts.back().empty()) break;
            
            // Add violated lazy contraints
            // TODO: unimplemented
            // for (; processedForts < forts.size(); ++processedForts) {
            //     GRBLinExpr fortSum;
            //     size_t blank = 0;
            //     for (auto& v: forts[processedForts]) {
            //         if (!state.graph().hasVertex(v)) {
            //             fmt::print("!!!invalid vertex {} in fort: {}!!!\n", v, forts[processedForts]);
            //         }
            //         if (state.isActive(v)) {
            //             if (output) {
            //                 fmt::print("???active vertex in neighborhood??? {} {} {}\n",
            //                            v,
            //                            lastSolution.numObserved(),
            //                            lastSolution.graph().numVertices());
            //             }
            //         }
            //         if (state.isBlank(v)) {
            //             fortSum += pi.at(v);
            //             ++blank;
            //         }
            //     }
            //     if (blank == 0) {
            //         fmt::print("??? infeasible fort: {} has no blank vertices\n", forts[processedForts]);
            //     }
            //     model.addConstr(fortSum >= 1);
            //     totalFortSize += forts[processedForts].size();
            //     if (output > 1) {
            //         fmt::print("fort {:4}: {} #{}({})\n", processedForts, forts[processedForts], forts[processedForts].size(), blank);
            //     }
            // }

            // Log cycle statics
            if (output) {
                fmt::print("#cycles: {}, avg size: {:.2f}\n", cycles.size(), double(totalCycleSize) / double(cycles.size()));
            }

            // Current time
            currentTime = now();

            // Remaining time for timeout
            remainingTimeout = std::max(0.0, timeLimit - std::chrono::duration_cast<std::chrono::duration<double>>(currentTime - startingTime).count());
            
            if (greedyUpper) {
                // Apply greedy heuristic to find a power dominating set
                fastGreedy(lastSolution, (greedyUpper - 1) * 2);
                size_t initial = lastSolution.numActive();
                topDownGreedy(lastSolution, false, blankVertices);
                if (lastSolution.numActive() < upperBound) {
                    upperBound = lastSolution.numActive();
                    if (output) {
                        fmt::print("greedy {} → {}\n", initial, lastSolution.numActive());
                    }
                    // Update feasible solution with the heuristic solution
                    std::swap(feasibleSolution, lastSolution);
                }
            }
            
            // Set time limit to the remaining time for timeout
            model.set(GRB_DoubleParam_TimeLimit, remainingTimeout);

            // Reoptimization
            model.optimize();

            // Get the status exit
            status = model.get(GRB_IntAttr_Status);
            
            // Objetive value of the last solution after reoptimization
            // It is a lower bound for the power domination number
            size_t new_bound = 0;

            switch (status) {

                // case GRB_INFEASIBLE: // Infeasible
                    // TODO: es necesario?
                    // return {1, 0, SolveState::Infeasible};

                // case GRB_INTERRUPTED: // Interrupted by the user
                    // TODO: es necesario?
                    // if (static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound)) < state.graph().numVertices()) {
                    //     new_bound = std::max(lowerBound, static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound)));
                    // }
                    // break;

                // case GRB_USER_OBJ_LIMIT: // User specified objective limit

                case GRB_OPTIMAL: // Optimal
                    // The new bound is the objetive value of the new (optimal) solution
                    new_bound = static_cast<size_t>(model.get(GRB_DoubleAttr_ObjVal) + 0.05);
                    break;

                case GRB_TIME_LIMIT: // Time limit
                    // The new bound is the final lower bound reported by gurobi
                    new_bound = static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound));
                    break;

                default:
                    fmt::print(stderr, "unexpected status: {}\n", status);
                    break;
            }

            if (model.get(GRB_IntAttr_SolCount) > 0) {
                // There is at least one solution

                //  Recover the new solution and update lastSolution
                for (auto v: state.graph().vertices()) {
                    if (sv.at(v).get(GRB_DoubleAttr_X) > 0.5) {
                        // If v is selected, then mark v as active 
                        lastSolution.setActive(v);
                    } else if (state.isBlank(v)) {
                        // if v is not selected and v is undecided in the previous solution,
                        // then mark v as undecided
                        lastSolution.setBlank(v);
                    } else {
                        // If is not selected and v is active in the previous solution,
                        // then mark v as inactive
                        lastSolution.setInactive(v);
                    }
                }
            }

            if (lowerBound > new_bound) {
                // Some checks
                if (status != GRB_TIME_LIMIT) {
                    fmt::print("!!!lowerBound decreased {} → {}!!!\n", lowerBound, new_bound);
                    if (size_t(model.get(GRB_DoubleAttr_ObjVal))) {
                        fmt::print("!!!wrong lower bound: obj {} < bound {}!!!\n", size_t(GRB_DoubleAttr_ObjVal), lowerBound);
                    }
                }
            } else if (new_bound <= state.graph().numVertices()) { // Only use valid bounds
                // Update the lower bound 
                lowerBound = new_bound;
            }

            if (lowerBound == upperBound) {
                // The new solution is an optimal power dominating set
                // Update feasible solution with the new solution
                std::swap(lastSolution, feasibleSolution);
                if (output) {
                    fmt::print("greedy solution is optimal\n");
                }
            } else if (lastSolution.allObserved()) {
                // The new solution is a feasible power dominating set
                // Update the upper bound
                // TODO: es necesario o ya lo hace el callback?
                upperBound = lastSolution.numActive();
            }

            // TODO: que hace???
            boundCallback(lowerBound, upperBound, cycles.size());

            if(earlyStop > 0) {
                // If earlyStop = 1, this stops the next reoptimization as soon as
                // a solution with objective value equal to lowerBound is found
                // (without finishing to prove of optimality).
                model.set(GRB_DoubleParam_BestObjStop, double(lowerBound));
            }
            
            if(output) { // Log statics
                fmt::print("LB: {}, UB: {} (status {}) (local LB: {}, UP: {})\n", lowerBound, upperBound, status, model.get(GRB_DoubleAttr_ObjBound), model.get(GRB_DoubleAttr_ObjVal));
            }
            
            if (remainingTimeout <= 1.0) status = GRB_TIME_LIMIT;
            
            if (lastSolution.allObserved()) {
                // The new solution is an optimal power dominating set
                // Update feasible solution with the new solution
                feasibleSolution = lastSolution;
                cycleCB(callback::When::FINAL, state, cycles, lowerBound, upperBound);
                break; // Finish loop
            } else {
                cycleCB(callback::When::INTERMEDIATE_HS, state, cycles, lowerBound, upperBound);
            }
            
            // For the next reoptimization, give gurobi the best power dominating set 
            // found so far as an initial solution
            // TODO: gurobi infiere valores para las demas variables?
            for (auto v: feasibleSolution.graph().vertices()) {
                if (feasibleSolution.isActive(v)) {
                    sv.at(v).set(GRB_DoubleAttr_Start, 1.0);
                } else {
                    sv.at(v).set(GRB_DoubleAttr_Start, 0.0);
                }
            }
        }

        // Update state with the optimal solution
        std::swap(state, feasibleSolution);
        if (upperBound == lowerBound) {
            status = GRB_OPTIMAL;
        }
        
        // TODO: ??
        if (upperBound < state.numActive()) {
            fmt::print("!!!upper bound too low!!!\n");
        }

        switch (status) {
            case GRB_OPTIMAL:
                return {lowerBound, lowerBound, SolveState::Optimal};
            case GRB_TIME_LIMIT:
                return {lowerBound, upperBound, SolveState::Timeout};
            default:
                return {size_t{0}, state.numActive() + state.numBlank(), SolveState::Timeout};
        }

    } catch (const GRBException& ex) {
        fmt::print(stderr, "Gurobi Error [{}]: {}\n", ex.getErrorCode(), ex.getMessage());
        throw ex;
    }
}

} //namespace pds
