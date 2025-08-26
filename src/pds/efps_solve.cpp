#include "efps_solve.hpp"

#include <cstdlib>
#include <ctime>
#include <gurobi_c++.h>
#include <range/v3/all.hpp>
#include <set>
#include <utility>

#include "gurobi_common.hpp"
#include "pdssolve.hpp"

namespace pds {

namespace {

// Function that gets the nearest visited predecessor of v, if any
std::pair<Node, bool> get_nearest_visited_predecessor(
    const PrecedenceDigraph &digraph,
    const std::map<Node, std::pair<Node, size_t>> &preceded_by, Node v) {
  Node pred;
  size_t max_level = 0;
  for (auto u : digraph.inNeighbors(v)) {
    if (!preceded_by.contains(u))
      continue;
    auto p = preceded_by.at(u);
    if (p.second <= max_level) {
      continue;
    }
    pred = u;
    max_level = p.second;
  }
  return std::make_pair(pred, max_level != 0);
}

// Function that gets the farthest visited successor of v different from u,
// if any
std::pair<Node, bool> get_farthest_visited_successor(
    const PrecedenceDigraph &digraph,
    const std::map<Node, std::pair<Node, size_t>> &preceded_by, Node v,
    Node u) {
  Node succ;
  size_t min_level = std::numeric_limits<size_t>::max();
  for (auto w : digraph.neighbors(v)) {
    if (!preceded_by.contains(w) || w == u)
      continue;
    auto p = preceded_by.at(w);
    if (p.second >= min_level) {
      continue;
    }
    succ = w;
    min_level = p.second;
  }
  return std::make_pair(succ, min_level != std::numeric_limits<size_t>::max());
}

// Functions that cuts a cycle from vertex u to vertex v
void cut_cycle(std::map<Node, std::pair<Node, size_t>> &preceded_by, Node u,
               Node v) {

  Node remove = u;
  Node pred;

  while (remove != v) {
    pred = preceded_by[remove].first;
    preceded_by.erase(remove);
    remove = pred;
  }

  return;
}

// Functions that cuts a cycle from vertex u to vertex v
// and makes u preceded by v
void cut_and_join_cycle(std::map<Node, std::pair<Node, size_t>> &preceded_by,
                        Node u, Node v) {

  size_t level = preceded_by[u].second;
  cut_cycle(preceded_by, u, v);
  preceded_by[u] = std::make_pair(v, level);

  return;
}

// Function that takes a cycle and makes it chordless
void make_cycle_chordless(const PrecedenceDigraph &digraph,
                          std::map<Node, std::pair<Node, size_t>> &preceded_by,
                          Node v) {
  Node pred = preceded_by[v].first;
  Node succ = v;
  // Check if pred has already visited successors (different from succ).
  // Select the farthest successor, i.e. with the lowest level.
  while (pred != v) {
    auto [u, ok] =
        get_farthest_visited_successor(digraph, preceded_by, pred, succ);
    if (ok)
      cut_and_join_cycle(preceded_by, u, pred);
    succ = pred;
    pred = preceded_by[succ].first;
  }
}

// Function that finds a chordless cycle in the precedence digraph, starting
// from vertex v (v is not guaranteed to be in the cycle).  We already know
// that such a cycle exists and that every vertex has at least one
// predecessor.
VertexList find_chordless_cycle(const PrecedenceDigraph &digraph, Node v) {
  // Map from vertex to predecessor and level (height)
  std::map<Node, std::pair<Node, size_t>> preceded_by;
  Node lastVertex = v;

  for (int i = 1; !preceded_by.contains(lastVertex); ++i) {
    assert(digraph.inDegree(lastVertex) > 0);

    // Check if lastVertex has already visited predecessors.
    // Select the nearest predecessor, i.e. with the highest level.
    // This strategy avoid chords in one direction, while the cycle must be
    // postprocessed to avoid chords in the other direction.
    auto [next, ok] =
        get_nearest_visited_predecessor(digraph, preceded_by, lastVertex);
    if (!ok)
      // Choose a random in-neighbor
      next = *(std::next(digraph.inNeighbors(lastVertex).begin(),
                         rand() % digraph.inDegree(lastVertex)));

    preceded_by.emplace(lastVertex, std::make_pair(next, i));
    lastVertex = next;
  }

  // Cut every vertex but the cycle
  cut_cycle(preceded_by, v, lastVertex);

  // Cut the cycle until getting a chordless cycle
  make_cycle_chordless(digraph, preceded_by, lastVertex);

  // Get cycle in terms of the original vertices
  VertexList cycle;
  auto u = lastVertex;
  do {
    cycle.push_back(u);
    u = preceded_by.at(u).first;
  } while (u != lastVertex);

  // Rotate the fps so the minium element is in the front
  ranges::rotate(cycle, ranges::min_element(cycle));

  return cycle;
}

// Function thats maps a cycle to an EFPS
EdgeList get_efps(VertexList &cycle, const EdgeMap<EdgeList> &prec2props) {
  EdgeList efps;
  for (auto it = cycle.rbegin(); it != cycle.rend();) {
    Node v = *it++;
    int u = it != cycle.rend() ? *it : cycle.back();
    for (auto &e : prec2props.at(v).at(u))
      efps.push_back(e);
  }
  return efps;
}

// Function that finds a set of EFPSs (associated to chordless cycles in the
// precedence digraph of the current solution)
// The size of the cycle is saved in the 2nd components.
std::set<std::pair<EdgeList, size_t>>
find_efpss(const PrecedenceDigraph &digraph,
           const EdgeMap<EdgeList> &prec2props, size_t fpssLimit) {
  std::set<std::pair<EdgeList, size_t>> efpss;
  // Iterate over the vertices
  for (auto v : digraph.vertices()) {
    // Check if the maximum number of FPSs has been found
    if (efpss.size() >= fpssLimit)
      break;
    // Find a chordless cycle (its existence is already guaranteed)
    VertexList cycle = find_chordless_cycle(digraph, v);
    // Transform the cycle into an EFPS
    EdgeList efps = get_efps(cycle, prec2props);
    efpss.insert(std::make_pair(efps, cycle.size()));
  }
  return efpss;
}

struct Callback : public GRBCallback {

  // Lower bound
  // TODO: que es
  size_t *lower;

  // An upper bound for the power dominating number
  // Number of active vertices + undecided vertices
  size_t *upper;

  // Map to GRBVar
  VertexMap<GRBVar> *sv;
  EdgeMap<GRBVar> *ye;

  // Base integer solution (before re-optimization)
  // It is not a power dominating set
  const PdsState *base;

  // Last integer solution (after re-optimizing)
  // It may or may not be a power dominating set
  PdsState *solution;

  // Last power dominating set found by gurobi
  PdsState *upperBound;

  // Span of undecided vertices
  std::span<PdsState::Vertex> blank;

  // TODO: ??
  int earlyStop;

  // Whether to add intermediate EFPS
  int intermediateEfps;

  // Vector of EFPS
  std::vector<std::pair<EdgeList, size_t>> *efps;

  // ???
  VertexSet *seen;

  // Verbosity
  int output;

  Callback(size_t &lower, size_t &upper, VertexMap<GRBVar> &sv,
           EdgeMap<GRBVar> &ye, const PdsState &base, PdsState &solution,
           PdsState &upperBound, std::span<PdsState::Vertex> blank,
           int earlyStop, int intermediateEfps,
           std::vector<std::pair<EdgeList, size_t>> &efps, VertexSet &seen,
           int output)
      : lower(&lower), upper(&upper), sv(&sv), ye(&ye), base(&base),
        solution(&solution), upperBound(&upperBound), blank(blank),
        earlyStop(earlyStop), intermediateEfps(intermediateEfps), efps(&efps),
        seen(&seen), output(output) {}

  void callback() override {

    // General MIP callback
    // Inside the NoRel Heuristic (when solving the root node and inside every
    // node)
    if (where == GRB_CB_MIP) {

      if (getIntInfo(GRB_CB_MIP_SOLCNT) > 0) {
        // At least one integer solution was found

        // MIP upper bound
        // TODO: GRB_CB_MIP_OBJBST or GRB_CB_MIPSOL_OBJ
        auto objVal =
            static_cast<size_t>(getDoubleInfo(GRB_CB_MIP_OBJBST) + 0.5);
        // MIP lower bound
        auto objBound = getDoubleInfo(GRB_CB_MIP_OBJBND);

        // If earlyStop = 2, this stops the reoptimization as soon as an invalid
        // solution is found, i.e. a non-power dominaintg set, whose objective
        // value is lower or equal to the upper bound (without finishing proving
        // optimality).
        // TODO: no entiendo: (objBound > *lower || earlyStop > 2)
        if (objVal <= *upper && !solution->allObserved() && earlyStop > 1 &&
            (objBound > *lower || earlyStop > 2)) {
          if (output) {
            fmt::print("early stop (exit 1) {} <= {}; {} <= {}\n", *lower,
                       objBound, objVal, *upper);
          }
          abort();
        }
      }
    }

    // MIP solution callback
    // An integer solution was found (it does not necessarily improve the
    // incumbent)
    if (where == GRB_CB_MIPSOL) {

      // Objective value of the new solution (rounded to the nearest integer)
      // static_cast truncates the number, so +0.5 is needed to round the number
      auto objVal = static_cast<size_t>(getDoubleInfo(GRB_CB_MIPSOL_OBJ) + 0.5);
      // MIP lower bound
      auto objBound = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);

      if (objVal <= *upper) {
        // The objective value of the solution is as good as the upper bound for
        // the power dominating number

        // Process the solution
        for (auto v : blank) {
          if (getSolution(sv->at(v)) > 0.5) {
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

        if (!solution->allObserved()) {
          // The solution is not a power dominating set

          // // TODO: parametros de violatedCycles ???
          // if (intermediateCycles >= 0) {
          //     // Cycles must be added for intermediate solutions
          //     for (auto& f: violatedCycles(*solution, intermediateCycles, 10,
          //     *seen, false)) {
          //         cycles->push_back(std::move(f));
          //     }
          // }

          if (output > 1) {
            fmt::print("new solution {} <= {}; {} <= {}; unobserved #{}\n",
                       *lower, objBound, objVal, *upper,
                       solution->numUnobserved());
          }

          // If earlyStop = 2, this stops the reoptimization as soon as an
          // invalid solution is found, i.e. a non-power dominaintg set, whose
          // objective value is lower or equal to the upper bound (without
          // finishing proving optimality).
          // TODO: por que se pide objBound > *lower
          if (earlyStop > 1 && objBound > *lower) {
            if (output) {
              fmt::print("early stop (exit 2) {} <= {}; {} <= {}\n", *lower,
                         objBound, objVal, *upper);
            }
            abort();
          }

        } else if (objVal < *upper) {
          // The solution is a power dominating set,
          // whose objetive value is better than the upper bound for the power
          // dominating number i.e. Gurobi improved the upper bound for the
          // power dominating number
          *upper = objVal;
          *upperBound = *solution;
          fmt::print("* new power dominating set #{}\n", objVal);
        }
      }
    }
  }
};

auto now() { return std::chrono::high_resolution_clock::now(); }

} // namespace

SolveResult solveEFPS(
    PdsState &state,  // Solution of the last reoptimization
    int output,       // Wheter to log
    double timeLimit, // Time limit
    int variant,      // TODO: ??
    int efpsInit,     // Wheter to initialize EFPS
    int greedyUpper,  // Whether to apply a greedy heuristic to find a pds
    int earlyStop,    // Whether to stop early
    efpscallback::EFPSCallback efpsCB, // Callback
    BoundCallback boundCallback,       // TODO: ??
    int intermediateEFPS, // Whether to add EFPS for intermediate solutions
    int efpsLimit         // Maximum number of EFPS add in each reoptimization
) {

  // Initialize lastSolution to be empty
  auto lastSolution = state;

  // Initialize feasibleSolution to be empty
  // Best power dominating set found so far
  // It can be updated by the MIP solution callback or the greedy heuristic
  PdsState feasibleSolution = state;

  // writePds(lastSolution.graph(), fmt::format("out/0_input.pds"));

  // Undecided vertices in the solution of the last reoptimization
  auto blankVertices =
      state.graph().vertices() |
      ranges::views::filter([&state](auto v) { return state.isBlank(v); }) |
      ranges::to<std::vector<PdsState::Vertex>>;

  // Set the lower bound to zero
  size_t lowerBound = 0;

  // The set of active or undecided vertices is a trivial power dominating set
  // This set gives an upper bound for the power dominating number
  // The upper bound can be updated by the MIP solution callback, the greedy
  // heuristic, or this function
  size_t upperBound = state.numActive() + state.numBlank();

  if (state.allObserved()) {
    // The solution of the last reoptimization is a power dominating set
    return {state.numActive(), state.numActive(), SolveState::Optimal};
  } else if (state.numBlank() == 0) {
    // The solution of the last reoptimization is not a power dominating set
    // and there is no undecided vertices
    return {state.graph().numVertices(), 0, SolveState::Infeasible};
  }

  // Set a random seed for the random generator
  srand((unsigned)time(0));

  // Map from precedences to propagations
  EdgeMap<EdgeList> prec2props;

  try {

    // TODO: ??
    VertexSet seen;

    // Initialize cycles
    // TODO: ???
    // auto cycles = initializeCycles(state, cycleInit, seen);
    std::vector<std::pair<EdgeList, size_t>> efpss;

    // TODO: remover estas lineas
    if (output > 1)
      fmt::print("Ignorar; cycleInit {} {}\n", efpsInit, variant);

    // Intitialize model
    GRBModel model(getEnv());

    // Set parameters
    model.set(GRB_IntParam_LogToConsole, false);
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    model.set(GRB_IntParam_NumericFocus, 1);
    model.set(GRB_DoubleParam_MIPGap, 1e-6);
    model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);

    // Set file for logging
    if (output) {
      model.set(GRB_StringParam_LogFile, "gurobi.log");
    }

    // Add variables
    VertexMap<GRBVar> sv;
    EdgeMap<GRBVar> ye;

    for (auto v : state.graph().vertices()) {

      // Variable: s_v
      if (state.isActive(v) || state.isBlank(v)) {
        sv.emplace(
            v, model.addVar(0.0, 1.0, 1.0, GRB_BINARY, fmt::format("s_{}", v)));
      }

      // Variable: y_e
      if (!state.isZeroInjection(v)) {
        continue;
      }

      for (auto e : state.graph().outEdges(v)) {
        auto u = state.graph().target(e);
        ye[v].emplace(u, model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                      fmt::format("y_{}{}", v, u)));
      }
    }

    // Add constraints

    // Constraints (1)
    // s_v + sum_{u in N(v)} (s_u + y_uv) >= 1, for all v in V
    for (auto v : state.graph().vertices()) {
      if (state.isActive(v)) {
        model.addConstr(sv.at(v) == 1.0);
      } else {
        GRBLinExpr observers = 0;
        if (state.isBlank(v)) {
          observers += sv.at(v);
        }
        for (auto e : state.graph().inEdges(v)) {
          auto u = state.graph().source(e);
          if (!state.isInactive(u))
            observers += sv.at(u);
          if (state.isZeroInjection(u))
            observers += ye.at(u).at(v);
        }
        model.addConstr(observers >= 1);
      }
    }

    // Limitation of outgoing propagations
    // (2) sum_{u \in N(v)} y_{vu} <= 1 - s_v, \forall v \in V_Z
    for (auto v : state.graph().vertices()) {
      if (!state.isZeroInjection(v)) {
        continue;
      }
      GRBLinExpr observers = 0;
      for (auto e : state.graph().outEdges(v)) {
        auto u = state.graph().target(e);
        observers += ye.at(v).at(u);
      }
      model.addConstr(observers <= 1 - sv.at(v));
    }

    // Build map from precedences to list of propagations
    for (auto v : state.graph().vertices()) {
      if (!state.isZeroInjection(v))
        continue;
      for (auto e : state.graph().outEdges(v)) {
        Node u = state.graph().target(e);
        prec2props[v][u].push_back(e);
        for (auto w : state.graph().neighbors(v)) {
          if (w == u)
            continue;
          prec2props[w][u].push_back(e);
        }
      }
    }

    // Add callback
    Callback cb(lowerBound, upperBound, sv, ye, state, lastSolution,
                feasibleSolution, blankVertices, earlyStop, intermediateEFPS,
                efpss, seen, output);
    model.setCallback(&cb);

    // Start the clock
    auto startingTime = now();

    // Set the initial status of the model to zero (unsolved)
    int status = 0;

    // Initialize statics for cycles
    size_t processedEfps = 0;
    size_t totalEfpsSize = 0;

    while (true) {

      // Current time
      auto currentTime = now();

      // Remaining time for timeout
      double remainingTimeout = std::max(
          0.0,
          timeLimit - std::chrono::duration_cast<std::chrono::duration<double>>(
                          currentTime - startingTime)
                          .count());

      if (model.get(GRB_IntAttr_SolCount) > 0) {
        // Add violated lazy contraints before reoptimization

        // Build precedence digraph (among unobserved vertices).
        // ATTENTION: redundant propagations are ignored and only
        // the unmonitored vertices are considered
        PrecedenceDigraph digraph;
        for (auto v : lastSolution.graph().vertices()) {
          if (lastSolution.isObserved(v)) {
            continue;
          }
          digraph.getOrAddVertex(v);
          bool redudant = false;
          for (auto u : lastSolution.graph().neighbors(v)) {
            if (!lastSolution.isZeroInjection(u)) {
              continue;
            }
            if (ye.at(u).at(v).get(GRB_DoubleAttr_X) < 0.5) {
              continue;
            }
            if (redudant)
              break;
            redudant = true;
            if (!lastSolution.isObserved(u)) {
              digraph.getOrAddVertex(u);
              digraph.addEdge(u, v);
            }
            for (auto w : lastSolution.graph().neighbors(u)) {
              if (w == v || lastSolution.isObserved(w)) {
                continue;
              }
              digraph.getOrAddVertex(w);
              digraph.addEdge(w, v);
            }
          }
        }

        // Find EFPSs
        auto moreEfpss = find_efpss(digraph, prec2props, efpsLimit);

        // Add the EFPSs to the set of EFPSs
        // Recall that the callback may have encountered some EFPSs in
        // intermediate solutions
        for (auto f : moreEfpss) {
          efpss.emplace_back(std::move(f));
        }

        // If there is no cycle, the instance is infeasible
        if (efpss.empty()) {
          return {state.graph().numVertices(), 0, SolveState::Infeasible};
        }

        // Add violated lazy contraints
        for (; processedEfps < efpss.size(); ++processedEfps) {
          GRBLinExpr efpsSum;
          auto &efps = efpss[processedEfps];
          for (auto e : efps.first)
            efpsSum += ye.at(digraph.source(e)).at(digraph.target(e));
          model.addConstr(efpsSum <= efps.second - 1);
          totalEfpsSize += efps.first.size();
          if (output > 1) {
            fmt::print("efps {:4}: {} #{}\n", processedEfps, efps.first,
                       efps.second);
          }
        }

        // Log cycle statics
        if (output) {
          fmt::print("#efps: {}, avg size: {:.2f}\n", efpss.size(),
                     double(totalEfpsSize) / double(efpss.size()));
        }
      }

      // Current time
      currentTime = now();

      // Remaining time for timeout
      remainingTimeout = std::max(
          0.0,
          timeLimit - std::chrono::duration_cast<std::chrono::duration<double>>(
                          currentTime - startingTime)
                          .count());

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
        // if (static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound)) <
        // state.graph().numVertices()) {
        //     new_bound = std::max(lowerBound,
        //     static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound)));
        // }
        // break;

        // case GRB_USER_OBJ_LIMIT: // User specified objective limit

      case GRB_OPTIMAL: // Optimal
        // The new bound is the objetive value of the new (optimal) solution
        new_bound = static_cast<size_t>(model.get(GRB_DoubleAttr_ObjVal) + 0.5);
        break;

      case GRB_TIME_LIMIT: // Time limit
        // The new bound is the final lower bound reported by gurobi
        new_bound = static_cast<size_t>(model.get(GRB_DoubleAttr_ObjBound));
        break;

      default:
        fmt::print(stderr, "unexpected status: {}\n", status);
        break;
      }

      if (model.get(GRB_IntAttr_SolCount) >
          0) { // There is at least one solution
        // lastSolution can differ from the last solution found in the MIPSOL
        // callback. Why? MIPSOL callback is only invoked when there is a new
        // incumbent solution (a solution with a better cost). So, the callback
        // won't be called for any other solution found after finding a solution
        // with optimal cost (gurobi might not know that the solution is optimal
        // until the dual bound is tightened). However, the optimal solutions
        // found later will be in the solution pool, and lastSolution will be
        // the last one.

        //  Recover the new solution and update lastSolution
        for (auto v : state.graph().vertices()) {
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
          fmt::print("!!!lowerBound decreased {} → {}!!!\n", lowerBound,
                     new_bound);
          if (size_t(model.get(GRB_DoubleAttr_ObjVal))) {
            fmt::print("!!!wrong lower bound: obj {} < bound {}!!!\n",
                       size_t(GRB_DoubleAttr_ObjVal), lowerBound);
          }
        }
      } else if (new_bound <=
                 state.graph().numVertices()) { // Only use valid bounds
        // Update the lower bound
        lowerBound = new_bound;
      }

      if (lowerBound == upperBound) {
        // The power dominating set in feasibleSolution is optimal (found in the
        // last MIPSOL callback) Update lastSolution to feasibleSolution
        std::swap(lastSolution, feasibleSolution);
        if (output) {
          fmt::print("power dominating set is optimal\n");
        }
      } else if (lastSolution.allObserved()) {
        // lastSolution is a power dominating set, so it is optimal
        // The upper bound can be updated
        fmt::print("* new power dominating set #{}\n",
                   lastSolution.numActive());
        fmt::print("power dominating set is optimal\n");
        upperBound = lastSolution.numActive();
      }

      // TODO: que hace???
      boundCallback(lowerBound, upperBound, efpss.size());

      if (earlyStop > 0) {
        // If earlyStop = 1, this stops the next reoptimization as soon as
        // a solution with objective value equal to lowerBound is found
        // (without finishing to prove of optimality).
        model.set(GRB_DoubleParam_BestObjStop, double(lowerBound));
      }

      if (output) { // Log statics
        fmt::print("LB: {}, UB: {} (status {}) (local LB: {}, UP: {})\n",
                   lowerBound, upperBound, status,
                   model.get(GRB_DoubleAttr_ObjBound),
                   model.get(GRB_DoubleAttr_ObjVal));
      }

      if (remainingTimeout <= 1.0)
        status = GRB_TIME_LIMIT;

      if (lastSolution.allObserved()) {
        // lastSolution is a power dominating set, so it is optimal
        // Update feasible solution with the new solution
        feasibleSolution = lastSolution;
        efpsCB(efpscallback::When::FINAL, state, efpss, lowerBound, upperBound);
        break; // Finish loop
      } else {
        efpsCB(efpscallback::When::INTERMEDIATE_HS, state, efpss, lowerBound,
               upperBound);
      }

      // For the next reoptimization, give gurobi the best power dominating set
      // found so far as an initial solution
      // TODO: gurobi infiere valores para las demas variables?
      for (auto v : feasibleSolution.graph().vertices()) {
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
      return {size_t{0}, state.numActive() + state.numBlank(),
              SolveState::Timeout};
    }

  } catch (const GRBException &ex) {
    fmt::print(stderr, "Gurobi Error [{}]: {}\n", ex.getErrorCode(),
               ex.getMessage());
    throw ex;
  }
}

} // namespace pds
