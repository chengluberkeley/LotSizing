//
//  LotSizingSolver.hpp
//  LotSizingSolver
//
//  Created by ChengLu on 11/24/21.
//  Copyright © 2021 Cheng Lu. All rights reserved.
//

#ifndef LotSizingSolver_hpp
#define LotSizingSolver_hpp

#include <vector>

/// Default: Linear cost edge
struct Edge {
    std::size_t from, to;
    int capacity;
    int flow;
    int cost;
};

// TODO: Convex cost edge

using ResidualEdge = Edge;

using productionEdges = std::vector<Edge>;

using forwardEdges = std::vector<Edge>;

using backwardEdges = std::vector<Edge>;

using residualEdges = std::vector<ResidualEdge>;

class ForwardGraph {
public:
    /// n: Number of demand nodes. Demand nodes are 1-indexed. The production node is always node 0. Thus the total number of nodes is n+1.
    ForwardGraph(std::size_t n, const std::vector<int>& productionCapacities, const std::vector<int>& productionCosts,
                 const std::vector<int>& forwardCapacities, const std::vector<int>& forwardCosts) {
        assert(n >= 1);
        assert(productionCapacities.size() == productionCosts.size());
        assert(n == productionCapacities.size());
        assert(forwardCapacities.size() == forwardCosts.size());
        assert(forwardCapacities.size() + 1 == n);

        // Production edges
        m_productionEdges.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            Edge& edge = m_productionEdges[i];
            edge.from = 0;
            edge.to = i + 1;
            edge.capacity = productionCapacities[i];
            edge.flow = 0;
            edge.cost = productionCosts[i];
        }

        // Inventory forward edges
        m_forwardEdges.resize(n - 1);
        for (std::size_t i = 0; i < n - 1; ++i) {
            Edge& edge = m_forwardEdges[i];
            edge.from = i + 1;
            edge.to = i + 2;
            edge.capacity = forwardCapacities[i];
            edge.flow = 0;
            edge.cost = forwardCosts[i];
        }
    }

protected:
    productionEdges m_productionEdges;
    forwardEdges m_forwardEdges;
};

class ForwardBackwardGraph: ForwardGraph {
public:
    ForwardBackwardGraph(std::size_t n, const std::vector<int>& productionCapacities, const std::vector<int>& productionCosts,
                         const std::vector<int>& forwardCapacities, const std::vector<int>& forwardCosts,
                         const std::vector<int>& backwardCapacities, const std::vector<int>& backwardCosts)
    : ForwardGraph(n, productionCapacities, productionCosts, forwardCapacities, forwardCosts)
    {
        assert(backwardCapacities.size() == backwardCosts.size());
        assert(backwardCapacities.size() == forwardCapacities.size());

        // Inventory backward edges
        m_backwardEdges.resize(n - 1);
        for (std::size_t i = 0; i < n - 1; ++i) {
            Edge& edge = m_backwardEdges[i];
            edge.from = i + 2;
            edge.to = i + 1;
            edge.capacity = backwardCapacities[i];
            edge.flow = 0;
            edge.cost = backwardCosts[i];
        }
    }

private:
    backwardEdges m_backwardEdges;
};

#endif /* LotSizingSolver_hpp */
