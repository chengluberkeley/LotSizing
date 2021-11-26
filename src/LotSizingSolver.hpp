//
//  LotSizingSolver.hpp
//  LotSizingSolver
//
//  Created by ChengLu on 11/24/21.
//  Copyright © 2021 Cheng Lu. All rights reserved.
//

#ifndef LotSizingSolver_hpp
#define LotSizingSolver_hpp

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
    void solve();

    bool constraintsSatisfied() const;

    bool isOptimal() const;

protected:
    std::size_t m_n;
    Demands m_demands;
    ProductionEdges m_productionEdges;
    InventoryEdges m_forwardEdges;

    // Forward residual edges
    ResidualEdges m_productionResidualEdges;
    ResidualEdges m_forwardResidualEdges;

    // Shared functions
    void orderedInsert(const ResidualPath& residualPath, std::list<ResidualPath>& residualPaths);

private:
    void elongateAndAdd(std::size_t node, std::list<ResidualPath>& residualPaths);

    void augmentAndUpdate(std::size_t node, std::list<ResidualPath>& residualPaths, uint32_t& demand);
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
    void solve();

private:
    InventoryEdges m_backwardEdges;

    // Backward residual edges
    ResidualEdges m_backwardResidualEdges;
};

#endif /* LotSizingSolver_hpp */
