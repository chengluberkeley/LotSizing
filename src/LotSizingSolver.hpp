//
//  LotSizingSolver.hpp
//  LotSizingSolver
//
//  Created by ChengLu on 11/24/21.
//  Copyright © 2021 Cheng Lu. All rights reserved.
//

#ifndef LotSizingSolver_hpp
#define LotSizingSolver_hpp

// #define DEBUG_LOTSIZING

#include "dp_array.h"

#include <list>
#include <vector>

/// Default: Linear cost edge
struct Edge {
    uint32_t capacity;
    uint32_t flow;
    double cost;
};

// TODO: Convex cost edge

using ProductionEdges = std::vector<Edge>;

using InventoryEdges = ProductionEdges;

using Demands = std::vector<uint32_t>;

// MARK: - Internal data structure used in the solver

using ResidualEdge = Edge;

using ResidualEdges = std::vector<ResidualEdge>;

/// Data for a residual path used in the O(n^2) algorithm.
struct ResidualPath {
    /// Starting node of the residual path.
    std::size_t from;

    /// Residual cost of the residual path.
    double cost;
};

// (start, end), both are inclusive.
using SegmentRange = std::pair<std::size_t, std::size_t>;

// (Non-zero capacitated backward residual paths ordered in costs, range of the residual paths)
using BackwardResidualPathSegment = std::pair<std::list<ResidualPath>, SegmentRange>;

class ForwardGraph {
public:
    /// n: Number of demand nodes. Demand nodes are 0-indexed. The production node is node n. Thus the total number of nodes is n+1.
    ForwardGraph(std::size_t n, const Demands& demands,
                 const std::vector<uint32_t>& productionCapacities,
                 const std::vector<double>& productionCosts,
                 const std::vector<uint32_t>& forwardCapacities,
                 const std::vector<double>& forwardCosts);

    /// Evaluate the total objective cost of the current solution at call time (may not be optimal).
    double cost() const;

    /// Solve the optimal solution to the problem instance.
    /// \return true if the problem is feasible.
    bool solve();

    /// Solve the optimal solution to the problem instance with faster algorithm.
    /// \return true if the problem is feasible.
    bool fasterSolve();

    /// Check whether the current solution satisfies all constraints.
    bool constraintsSatisfied() const;

    /// Check whether the current solution is optimal.
    bool isOptimal() const;

protected:
    std::size_t m_n;
    Demands m_demands;
    ProductionEdges m_productionEdges;
    InventoryEdges m_forwardEdges;
    double m_costSum = 0;
    double m_costBound = 0;

    // Production and forward residual edges
    ResidualEdges m_productionResidualEdges;
    ResidualEdges m_forwardResidualEdges;

    // Shared functions
    void orderedInsert(const ResidualPath& residualPath, std::list<ResidualPath>& residualPaths);

    bool isInfinite(double value);

#ifdef DEBUG_LOTSIZING
    void print() const;

    void printResiduals() const;
#endif

private:
    void elongateAndAdd(std::size_t node, std::list<ResidualPath>& residualPaths);

    // Return true if the augmentation is successful.
    bool augmentAndUpdate(std::size_t node, std::list<ResidualPath>& residualPaths, uint32_t& demand);

    // MARK: - Dynamic path version

    void fastElongateAndAdd(std::size_t node, std::size_t& start,
                            dp_array<uint32_t>& capacityDP, dp_array<double>& costDP);

    // Return true if the augmentation is successful.
    bool fastAugmentAndUpdate(std::size_t node, std::size_t& start,
                              dp_array<uint32_t>& capacityDP, dp_array<double>& costDP, dp_array<uint32_t>& flowDP,
                              uint32_t& demand);
};

class ForwardBackwardGraph: ForwardGraph {
public:
    ForwardBackwardGraph(std::size_t n, const Demands& demands,
                         const std::vector<uint32_t>& productionCapacities,
                         const std::vector<double>& productionCosts,
                         const std::vector<uint32_t>& forwardCapacities,
                         const std::vector<double>& forwardCosts,
                         const std::vector<uint32_t>& backwardCapacities,
                         const std::vector<double>& backwardCosts);

    /// Evaluate the total objective cost of the current solution at call time (may not be optimal).
    double cost() const;

    /// Solve the optimal solution to the problem instance.
    /// \return true if the problem is feasible.
    bool solve();

    /// Solve the optimal solution to the problem instance with faster algorithm.
    /// \return true if the problem is feasible.
    bool fasterSolve();

    /// Check whether the current solution satisfies all constraints.
    bool constraintsSatisfied() const;

    /// Check whether the current solution is optimal.
    bool isOptimal() const;

private:
    InventoryEdges m_backwardEdges;

    // Backward residual edges
    ResidualEdges m_backwardResidualEdges;

    // Auxiliary data structure
    std::vector<double> m_backwardResidualEdgeCostSums;

    void elongateAndUpdate(std::size_t node, std::list<ResidualPath>& forwardResidualPaths,
                           std::list<BackwardResidualPathSegment>& backwardResidualPathSegments);

    // Return true if the augmentation is successful.
    bool augmentAndUpdate(std::size_t node, std::list<ResidualPath>& forwardResidualPaths,
                          std::list<BackwardResidualPathSegment>& backwardResidualPathSegments,
                          uint32_t& demand);

    // MARK: - Dynamic path version

    void fastElongateAndUpdate(std::size_t node, const dp_array<uint32_t>& forwardResFlowDP, const dp_array<uint32_t>& backwardResFlowDP,
                               std::size_t& start, std::list<SegmentRange>& backwardResidualSegments,
                               dp_array<uint32_t>& forwardResCapacityDP, dp_array<uint32_t>& backwardResCapacityDP,
                               dp_array<double>& costDP);

    // Return true if the augmentation is successful.
    bool fastAugmentAndUpdate(std::size_t node, std::size_t& start, std::list<SegmentRange>& backwardResidualSegments,
                              dp_array<uint32_t>& forwardResCapacityDP, dp_array<uint32_t>& backwardResCapacityDP,
                              dp_array<double>& costDP,
                              dp_array<uint32_t>& forwardResFlowDP, dp_array<uint32_t>& backwardResFlowDP,
                              uint32_t& demand);

    // Convert residual flows in DPs to final solution.
    void mergeFlowSolutions(const dp_array<uint32_t>& forwardResFlowDP, const dp_array<uint32_t>& backwardResFlowDP);

#ifdef DEBUG_LOTSIZING
    void print() const;

    void print(const std::list<ResidualPath>& forwardResidualPaths,
               const std::list<BackwardResidualPathSegment>& backwardResidualPathSegments) const;

    void printResiduals() const;
#endif
};

#endif /* LotSizingSolver_hpp */
