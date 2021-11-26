//
//  LotSizingSolver.cpp
//  LotSizingSolver
//
//  Created by ChengLu on 11/24/21.
//  Copyright © 2021 Cheng Lu. All rights reserved.
//

#include "LotSizingSolver.hpp"

// MARK: - ForwardGraph functions
ForwardGraph::ForwardGraph(std::size_t n, const Demands& demands,
                           const std::vector<uint32_t>& productionCapacities,
                           const std::vector<double>& productionCosts,
                           const std::vector<uint32_t>& forwardCapacities,
                           const std::vector<double>& forwardCosts) {
    assert(n >= 1);
    assert(demands.size() == n);
    assert(productionCapacities.size() == productionCosts.size());
    assert(productionCapacities.size() == n);
    assert(forwardCapacities.size() == forwardCosts.size());
    assert(forwardCapacities.size() + 1 == n);

    m_n = n;

    // Demands & Production edges
    uint32_t totalDemand = 0;
    m_demands.resize(n);
    uint32_t totalProductionCapacity = 0;
    m_productionEdges.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        m_demands[i] = demands[i];
        totalDemand += demands[i];

        Edge& edge = m_productionEdges[i];
        edge.capacity = productionCapacities[i];
        edge.flow = 0;
        edge.cost = productionCosts[i];
        assert(edge.cost >= 0);
        totalProductionCapacity += edge.capacity;
        // Partial-sum satisfactory as well.
        assert(totalDemand <= totalProductionCapacity);
    }

    // Inventory forward edges
    m_forwardEdges.resize(n - 1);
    for (std::size_t i = 0; i < n - 1; ++i) {
        Edge& edge = m_forwardEdges[i];
        edge.capacity = forwardCapacities[i];
        edge.flow = 0;
        edge.cost = forwardCosts[i];
        assert(edge.cost >= 0);
    }

    m_productionResidualEdges = m_productionEdges;
    m_forwardResidualEdges = m_forwardEdges;
}

double ForwardGraph::cost() const {
    double cost = 0;
    for (const auto& edge: m_productionEdges) {
        cost += edge.cost * edge.flow;
    }

    for (const auto& edge: m_forwardEdges) {
        cost += edge.cost * edge.flow;
    }

    return cost;
}

void ForwardGraph::solve() {
    std::list<ResidualPath> residualPaths;

    for (std::size_t i = 0; i < m_n; ++i) {
        elongateAndAdd(i, residualPaths);
        uint32_t demand = m_demands[i];
        while (demand > 0) {
            augmentAndUpdate(i, residualPaths, demand);
        }
    }
}

bool ForwardGraph::constraintsSatisfied() const {
    // Capacity constraints
    for (const auto& edge: m_productionEdges) {
        if (edge.flow > edge.capacity) {
            return false;
        }
    }

    for (const auto& edge: m_forwardEdges) {
        if (edge.flow > edge.capacity) {
            return false;
        }
    }

    // Flow balance constraints
    // First node
    if (m_productionEdges[0].flow - m_forwardEdges[0].flow != m_demands[0]) {
        return false;
    }

    for (std::size_t i = 1; i < m_n - 1; ++i) {
        if (m_productionEdges[i].flow + m_forwardEdges[i - 1].flow - m_forwardEdges[i].flow != m_demands[i]) {
            return false;
        }
    }

    if (m_productionEdges[m_n - 1].flow + m_forwardEdges[m_n - 2].flow != m_demands[m_n - 1]) {
        return false;
    }

    return true;
}

bool ForwardGraph::isOptimal() const {
    std::vector<double> d(m_n, std::numeric_limits<double>::max());

    // Bellman-Ford shortest path algorithm.

    for (std::size_t i = 0; i < m_n - 1; ++i) {
        // Production arcs
        for (std::size_t j = 0; j < m_n; ++j) {
            const auto& edge = m_productionEdges[j];
            if (edge.flow < edge.capacity) {
                if (edge.cost < d[j]) {
                    d[j] = edge.cost;
                }
            }

            if (edge.flow > 0) {
                if (d[j] - edge.cost < 0) {
                    return false;
                }
            }
        }

        for (std::size_t j = 0; j < m_n - 1; ++j) {
            const auto& edge = m_forwardEdges[j];
            if (edge.flow < edge.capacity) {
                if (d[j] + edge.cost < d[j + 1]) {
                    d[j + 1] = d[j] + edge.cost;
                }
            }

            if (edge.flow > 0) {
                if (d[j + 1] - edge.cost < d[j]) {
                    d[j] = d[j + 1] - edge.cost;
                }
            }
        }
    }

    // Verify no-negative cycles.
    for (std::size_t i = 0; i < m_n; ++i) {
        const auto& edge = m_productionEdges[i];
        if (edge.flow < edge.capacity) {
            if (edge.cost < d[i]) {
                return false;
            }
        }

        if (edge.flow > 0) {
            if (d[i] - edge.cost < 0) {
                return false;
            }
        }
    }

    for (std::size_t i = 0; i < m_n - 1; ++i) {
        const auto& edge = m_forwardEdges[i];
        if (edge.flow < edge.capacity) {
            if (d[i] + edge.cost < d[i + 1]) {
                return false;
            }
        }

        if (edge.flow > 0) {
            if (d[i + 1] - edge.cost < d[i]) {
                return false;
            }
        }
    }

    return true;
}

// MARK: - Protected

void ForwardGraph::orderedInsert(const ResidualPath& residualPath, std::list<ResidualPath>& residualPaths) {
    auto it = residualPaths.begin();
    while (it != residualPaths.end() && residualPath.cost > it->cost) {
        ++it;
    }

    if (it == residualPaths.end()) {
        residualPaths.push_back(residualPath);
        return;
    }

    residualPaths.insert(it, residualPath);
}

// MARK: - Private

void ForwardGraph::elongateAndAdd(std::size_t node, std::list<ResidualPath>& residualPaths) {
    if (node > 0) {
        // First append arc (node - 1, node) to each existing residual path.
        const ResidualEdge& residualEdge = m_forwardResidualEdges[node - 1];

        if (residualEdge.capacity > 0) {
            for (auto it = residualPaths.begin(); it != residualPaths.end(); ++it) {
                it->cost += residualEdge.cost;
            }

            // Insert the residual path (0) -> (node-1) -> (node)
            if (m_productionResidualEdges[node - 1].capacity > 0) {
                ResidualPath residualPath;
                residualPath.from = node - 1;
                residualPath.cost = m_productionResidualEdges[node - 1].cost + residualEdge.cost;
                // Insert
                orderedInsert(residualPath, residualPaths);
            }
        } else {
            // Clear all existing residual paths.
            residualPaths.clear();
        }
    }

    // Insert the residual path (0) -> (node)
    if (m_productionResidualEdges[node].capacity > 0) {
        ResidualPath residualPath;
        residualPath.from = node;
        residualPath.cost = m_productionResidualEdges[node].cost;
        // Insert
        orderedInsert(residualPath, residualPaths);
    }
}

void ForwardGraph::augmentAndUpdate(std::size_t node, std::list<ResidualPath>& residualPaths, uint32_t& demand) {
    assert(demand > 0);

    // Get the minimum cost residual path.
    const auto& residualPath = residualPaths.front();

    // Find residual capacity
    uint32_t residualCapacity = m_productionResidualEdges[residualPath.from].capacity;
    for (std::size_t i = residualPath.from; i < node; ++i) {
        residualCapacity = std::min(residualCapacity, m_forwardResidualEdges[i].capacity);
    }
    assert(residualCapacity > 0);

    auto delta = std::min(demand, residualCapacity);
    assert(delta > 0);

    // Update demand
    demand -= delta;

    // Update flows and residual edges.
    // TODO: Re-visit to save unnecessary computation.
    std::size_t lastZeroResCapNode = m_n;
    assert(m_productionResidualEdges[residualPath.from].capacity >= delta);
    m_productionResidualEdges[residualPath.from].capacity -= delta;
    m_productionEdges[residualPath.from].flow += delta;
    for (std::size_t i = residualPath.from; i < node; ++i) {
        assert(m_forwardResidualEdges[i].capacity >= delta);
        m_forwardResidualEdges[i].capacity -= delta;
        // For forward graph, the forward edges are the same as the forward residual edges.
        m_forwardEdges[i].flow += delta;
        if (m_forwardResidualEdges[i].capacity == 0) {
            lastZeroResCapNode = i;
        }
    }

    // Update the residual paths.
    residualPaths.pop_front();

    if (residualPath.from == node) {
        return;
    }

    if (lastZeroResCapNode < m_n) {
        residualPaths.remove_if([&node, &lastZeroResCapNode](const ResidualPath& residualPath) {
            return residualPath.from != node && residualPath.from <= lastZeroResCapNode;
        });
    }
}

// MARK: - ForwardBackwardGraph functions

ForwardBackwardGraph::ForwardBackwardGraph(std::size_t n, const Demands& demands,
                                           const std::vector<uint32_t>& productionCapacities,
                                           const std::vector<double>& productionCosts,
                                           const std::vector<uint32_t>& forwardCapacities,
                                           const std::vector<double>& forwardCosts,
                                           const std::vector<uint32_t>& backwardCapacities,
                                           const std::vector<double>& backwardCosts)
: ForwardGraph(n, demands, productionCapacities, productionCosts, forwardCapacities, forwardCosts) {
    assert(backwardCapacities.size() == backwardCosts.size());
    assert(backwardCapacities.size() == forwardCapacities.size());

    // Inventory backward edges
    m_backwardEdges.resize(n - 1);
    for (std::size_t i = 0; i < n - 1; ++i) {
        Edge& edge = m_backwardEdges[i];
        edge.capacity = backwardCapacities[i];
        edge.flow = 0;
        edge.cost = backwardCosts[i];
        assert(edge.cost >= 0);
    }

    m_backwardResidualEdges = m_backwardEdges;
}

double ForwardBackwardGraph::cost() const {
    double cost = ForwardGraph::cost();

    for (const auto& edge: m_backwardEdges) {
        cost += edge.cost * edge.flow;
    }

    return cost;
}

void ForwardBackwardGraph::solve() {
    // TODO:
}
