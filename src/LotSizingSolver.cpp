//
//  LotSizingSolver.cpp
//  LotSizingSolver
//
//  Created by ChengLu on 11/24/21.
//  Copyright © 2021 Cheng Lu. All rights reserved.
//

#include "LotSizingSolver.hpp"

#include <algorithm>
#include <utility>

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
    m_backwardResidualEdgeCostSums.resize(n - 1);
    for (std::size_t i = 0; i < n - 1; ++i) {
        Edge& edge = m_backwardEdges[i];
        edge.capacity = backwardCapacities[i];
        edge.flow = 0;
        edge.cost = backwardCosts[i];
        assert(edge.cost >= 0);
        if (i == 0) {
            m_backwardResidualEdgeCostSums[i] = edge.cost;
        } else {
            m_backwardResidualEdgeCostSums[i] = m_backwardResidualEdgeCostSums[i - 1] + edge.cost;
        }
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
    std::list<ResidualPath> forwardResidualPaths;
    std::list<BackwardResidualPathSegment> backwardResidualPathSegments;

    // Initialize the forward and backward residual paths.
    if (m_productionResidualEdges[0].capacity > 0) {
        ResidualPath residualPath;
        residualPath.from = 0;
        residualPath.cost = m_productionResidualEdges[0].cost;
        orderedInsert(residualPath, forwardResidualPaths);
    }

    // Accumulated backward residual cost.
    std::list<ResidualPath> backwardResidualPath;
    std::size_t start = 0;
    std::size_t end = start;
    for (std::size_t i = 1; i < m_n; ++i) {
        if (m_backwardResidualEdges[i - 1].capacity == 0) {
            if (start < end && !backwardResidualPath.empty()) {
                backwardResidualPathSegments.push_back(std::make_pair(backwardResidualPath,
                                                                      std::make_pair(start, end)));
            }
            start = i;
            end = start;
            backwardResidualPath.clear();
        } else {
            // Ordered insert a backward residual path.
            if (m_productionResidualEdges[i].capacity > 0 && i != start) {
                ResidualPath residualPath;
                residualPath.from = i;
                residualPath.cost = m_backwardResidualEdgeCostSums[i - 1] + m_productionResidualEdges[i].cost;
                orderedInsert(residualPath, backwardResidualPath);
            }
            end = i;
        }
    }

    // Last segment
    assert(end == m_n - 1);
    if (start < end && !backwardResidualPath.empty()) {
        backwardResidualPathSegments.push_back(std::make_pair(backwardResidualPath,
                                                              std::make_pair(start, end)));
    }

    for (std::size_t i = 0; i < m_n; ++i) {
        uint32_t demand = m_demands[i];
        while (demand > 0) {
            augmentAndUpdate(i, forwardResidualPaths, backwardResidualPathSegments, demand);
        }
        if (i < m_n - 1) {
            elongateAndUpdate(i + 1, forwardResidualPaths, backwardResidualPathSegments);
        }
    }
}

bool ForwardBackwardGraph::constraintsSatisfied() const {
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

    for (const auto& edge: m_backwardEdges) {
        if (edge.flow > edge.capacity) {
            return false;
        }
    }

    for (std::size_t i = 0; i < m_n - 1; ++i) {
        if (m_forwardEdges[i].flow > 0 && m_backwardEdges[i].flow > 0) {
            return false;
        }
    }

    // Flow balance constraints
    // First node
    if (m_productionEdges[0].flow + m_backwardEdges[0].flow - m_forwardEdges[0].flow
        != m_demands[0]) {
        return false;
    }

    for (std::size_t i = 1; i < m_n - 1; ++i) {
        if (m_productionEdges[i].flow + m_forwardEdges[i - 1].flow + m_backwardEdges[i].flow
            - m_forwardEdges[i].flow - m_backwardEdges[i - 1].flow
            != m_demands[i]) {
            return false;
        }
    }

    if (m_productionEdges[m_n - 1].flow + m_forwardEdges[m_n - 2].flow - m_backwardEdges[m_n - 2].flow
        != m_demands[m_n - 1]) {
        return false;
    }

    return true;
}

bool ForwardBackwardGraph::isOptimal() const {
    std::vector<double> d(m_n, std::numeric_limits<double>::max());

    for (std::size_t j = 0; j < m_n - 1; ++j) {
        const auto& fEdge = m_forwardResidualEdges[j];
        if (fEdge.capacity > 0) {
            assert(!(m_forwardEdges[j].flow > 0 && m_backwardEdges[j].flow > 0));
            if (m_forwardEdges[j].flow == 0 && m_backwardEdges[j].flow == 0) {
                assert(fEdge.capacity == m_forwardEdges[j].capacity);
                assert(fEdge.cost == m_forwardEdges[j].cost);
            } else if (m_forwardEdges[j].flow > 0) {
                assert(fEdge.capacity == m_forwardEdges[j].capacity - m_forwardEdges[j].flow);
                assert(fEdge.cost == m_forwardEdges[j].cost);
            } else if (m_backwardEdges[j].flow > 0) {
                assert(fEdge.capacity == m_backwardEdges[j].flow);
                assert(fEdge.cost == -m_backwardEdges[j].cost);
            }
        }

        const auto& bEdge = m_backwardResidualEdges[j];
        if (bEdge.capacity > 0) {
            assert(!(m_forwardEdges[j].flow > 0 && m_backwardEdges[j].flow > 0));
            if (m_forwardEdges[j].flow == 0 && m_backwardEdges[j].flow == 0) {
                assert(bEdge.capacity == m_backwardEdges[j].capacity);
                assert(bEdge.cost == m_backwardEdges[j].cost);
            } else if (m_backwardEdges[j].flow > 0) {
                assert(bEdge.capacity == m_backwardEdges[j].capacity - m_backwardEdges[j].flow);
                assert(bEdge.cost == m_backwardEdges[j].cost);
            } else if (m_forwardEdges[j].flow > 0) {
                assert(bEdge.capacity == m_forwardEdges[j].flow);
                assert(bEdge.cost == -m_forwardEdges[j].cost);
            }
        }
    }

    // Bellman-Ford shortest path algorithm.

    for (std::size_t i = 0; i < m_n + 1; ++i) {
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
            const auto& fEdge = m_forwardResidualEdges[j];
            if (fEdge.capacity > 0) {
                assert(!(m_forwardEdges[j].flow > 0 && m_backwardEdges[j].flow > 0));
                if (m_forwardEdges[j].flow == 0 && m_backwardEdges[j].flow == 0) {
                    assert(fEdge.capacity == m_forwardEdges[j].capacity);
                    assert(fEdge.cost == m_forwardEdges[j].cost);
                } else if (m_forwardEdges[j].flow > 0) {
                    assert(fEdge.capacity == m_forwardEdges[j].capacity - m_forwardEdges[j].flow);
                    assert(fEdge.cost == m_forwardEdges[j].cost);
                } else if (m_backwardEdges[j].flow > 0) {
                    assert(fEdge.capacity == m_backwardEdges[j].flow);
                    assert(fEdge.cost == -m_backwardEdges[j].cost);
                }
                if (d[j] + fEdge.cost < d[j + 1]) {
                    d[j + 1] = d[j] + fEdge.cost;
                }
            }

            const auto& bEdge = m_backwardResidualEdges[j];
            if (bEdge.capacity > 0) {
                assert(!(m_forwardEdges[j].flow > 0 && m_backwardEdges[j].flow > 0));
                if (m_forwardEdges[j].flow == 0 && m_backwardEdges[j].flow == 0) {
                    assert(bEdge.capacity == m_backwardEdges[j].capacity);
                    assert(bEdge.cost == m_backwardEdges[j].cost);
                } else if (m_backwardEdges[j].flow > 0) {
                    assert(bEdge.capacity == m_backwardEdges[j].capacity - m_backwardEdges[j].flow);
                    assert(bEdge.cost == m_backwardEdges[j].cost);
                } else if (m_forwardEdges[j].flow > 0) {
                    assert(bEdge.capacity == m_forwardEdges[j].flow);
                    assert(bEdge.cost == -m_forwardEdges[j].cost);
                }
                if (d[j + 1] + bEdge.cost < d[j]) {
                    d[j] = d[j + 1] + bEdge.cost;
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
        const auto& fEdge = m_forwardResidualEdges[i];
        if (fEdge.capacity > 0) {
            if (d[i] + fEdge.cost < d[i + 1]) {
                return false;
            }
        }

        const auto& bEdge = m_backwardResidualEdges[i];
        if (bEdge.capacity > 0) {
            if (d[i + 1] + bEdge.cost < d[i]) {
                return false;
            }
        }
    }

    return true;
}

// MARK: - Private

void ForwardBackwardGraph::elongateAndUpdate(std::size_t node, std::list<ResidualPath>& forwardResidualPaths,
                                             std::list<BackwardResidualPathSegment>& backwardResidualPathSegments) {
    assert(node > 0);

    // MARK: - Forward residual paths.

    // First append residual arc (node - 1, node) to each existing forward residual path.
    const auto& residualEdge = m_forwardResidualEdges[node - 1];

    if (residualEdge.capacity > 0) {
        for (auto it = forwardResidualPaths.begin(); it != forwardResidualPaths.end(); ++it) {
            it->cost += residualEdge.cost;
        }
    } else {
        // Clear all existing forward residual paths.
        forwardResidualPaths.clear();
    }

    // Insert the residual path (0) -> (node)
    if (m_productionResidualEdges[node].capacity > 0) {
        ResidualPath residualPath;
        residualPath.from = node;
        residualPath.cost = m_productionResidualEdges[node].cost;
        // Insert
        orderedInsert(residualPath, forwardResidualPaths);
    }

    // MARK: - Backward residual paths.

    if (backwardResidualPathSegments.empty()) {
        return;
    }

    auto& backwardResidualPathSegment = backwardResidualPathSegments.front();
    assert(!backwardResidualPathSegment.first.empty());
    if (backwardResidualPathSegment.second.first > node) {
        return;
    }

    assert(backwardResidualPathSegment.second.first >= node - 1
           && backwardResidualPathSegment.second.first <= node);
    backwardResidualPathSegment.second.first = node;
    if (backwardResidualPathSegment.second.first >= backwardResidualPathSegment.second.second) {
        backwardResidualPathSegments.pop_front();
        return;
    }

    backwardResidualPathSegment.first.remove_if([&node](const ResidualPath& resPath) {
        return resPath.from <= node;
    });

    if (backwardResidualPathSegment.first.empty()) {
        backwardResidualPathSegments.pop_front();
    }
}

void ForwardBackwardGraph::augmentAndUpdate(std::size_t node, std::list<ResidualPath>& forwardResidualPaths,
                                            std::list<BackwardResidualPathSegment>& backwardResidualPathSegments,
                                            uint32_t& demand) {
    // TODO: If both are empty, it indicates non-feasibility.
    assert(!forwardResidualPaths.empty() || !backwardResidualPathSegments.empty());

    ResidualPath residualPath;

    // Get the minimum cost out of forward and backward residual paths.
    if (!forwardResidualPaths.empty()) {
        residualPath = forwardResidualPaths.front();
    }

    if (!backwardResidualPathSegments.empty()) {
        const auto& backwardResidualPathSegment = backwardResidualPathSegments.front();
        const auto firstNode = backwardResidualPathSegment.second.first;
        assert(node <= firstNode);
        if (node == firstNode) {
            assert(!backwardResidualPathSegment.first.empty());
            const auto& backwardResidualPath = backwardResidualPathSegment.first.front();

            double cost = backwardResidualPath.cost;
            if (firstNode > 0) {
                cost -= m_backwardResidualEdgeCostSums[firstNode - 1];
            }

            if (forwardResidualPaths.empty() || cost < residualPath.cost) {
                residualPath = backwardResidualPath;
            }
        }
    }

    // Find residual capacity
    uint32_t residualCapacity = m_productionResidualEdges[residualPath.from].capacity;

    if (node >= residualPath.from) {
        // Forward residual path
        for (std::size_t i = residualPath.from; i < node; ++i) {
            residualCapacity = std::min(residualCapacity, m_forwardResidualEdges[i].capacity);
        }
        assert(residualCapacity > 0);

        auto delta = std::min(demand, residualCapacity);
        assert(delta > 0);

        // Update demand
        demand -= delta;

        // Update flows and residual edges
        std::size_t lastZeroResCapNode = m_n;
        std::vector<std::pair<std::size_t, double>> residualCostUpdatedEdges;
        bool producationSaturated = false;

        assert(m_productionResidualEdges[residualPath.from].capacity >= delta);
        m_productionResidualEdges[residualPath.from].capacity -= delta;
        m_productionEdges[residualPath.from].flow += delta;
        assert(m_productionEdges[residualPath.from].flow <= m_productionEdges[residualPath.from].capacity);
        if (m_productionResidualEdges[residualPath.from].capacity == 0) {
            producationSaturated = true;
        }
        for (std::size_t i = residualPath.from; i < node; ++i) {
            assert(m_forwardResidualEdges[i].capacity >= delta);
            m_forwardResidualEdges[i].capacity -= delta;
            // Two cases
            if (m_backwardEdges[i].flow > 0) {
                // Reverse arc
                assert(m_backwardEdges[i].flow >= delta);
                assert(m_forwardEdges[i].flow == 0);
                m_backwardEdges[i].flow -= delta;

                m_backwardResidualEdges[i].capacity = m_backwardEdges[i].capacity - m_backwardEdges[i].flow;
                if (m_forwardResidualEdges[i].capacity == 0) {
                    assert(m_backwardEdges[i].flow == 0);
                    // Turn to the regular arc
                    if (m_forwardEdges[i].capacity > 0) {
                        double costDelta = m_backwardEdges[i].cost + m_forwardEdges[i].cost;
                        residualCostUpdatedEdges.push_back(std::make_pair(i, costDelta));
                        m_forwardResidualEdges[i].capacity = m_forwardEdges[i].capacity;
                        m_forwardResidualEdges[i].cost = m_forwardEdges[i].cost;
                    } else {
                        lastZeroResCapNode = i;
                    }
                }
            } else {
                // Regular arc
                m_forwardEdges[i].flow += delta;
                m_backwardResidualEdges[i].capacity = m_forwardEdges[i].flow;
                m_backwardResidualEdges[i].cost = -m_forwardEdges[i].cost;
                assert(m_forwardEdges[i].flow <= m_forwardEdges[i].capacity);
                if (m_forwardResidualEdges[i].capacity == 0) {
                    lastZeroResCapNode = i;
                }
            }
        }

        // Update the forward residual paths.
        if (producationSaturated || lastZeroResCapNode < m_n) {
            forwardResidualPaths.pop_front();
        }

        if (residualPath.from == node) {
            return;
        }

        if (lastZeroResCapNode < m_n) {
            forwardResidualPaths.remove_if([&node, &lastZeroResCapNode](const ResidualPath& residualPath) {
                assert(residualPath.from < node);
                return residualPath.from <= lastZeroResCapNode;
            });
        }

        if (!residualCostUpdatedEdges.empty() && !forwardResidualPaths.empty()) {
            std::for_each(forwardResidualPaths.begin(), forwardResidualPaths.end(),
                          [&residualCostUpdatedEdges](ResidualPath& residualPath) {
                int index = static_cast<int>(residualCostUpdatedEdges.size()) - 1;
                while (index >= 0 && residualPath.from <= residualCostUpdatedEdges[index].first) {
                    residualPath.cost += residualCostUpdatedEdges[index].second;
                    --index;
                }
            });
            forwardResidualPaths.sort([](const ResidualPath& resPath0, const ResidualPath& resPath1) {
                return resPath0.cost < resPath1.cost;
            });
        }
    } else {
        assert(residualPath.from > node);

        // Backward residual path
        for (std::size_t i = node; i < residualPath.from; ++i) {
            residualCapacity = std::min(residualCapacity, m_backwardResidualEdges[i].capacity);
        }
        assert(residualCapacity > 0);

        auto delta = std::min(demand, residualCapacity);
        assert(delta > 0);

        // Update demand
        demand -= delta;

        // Update flows and residual edges
        std::vector<std::size_t> zeroBackwardResidualEdges;
        bool productionSaturated = false;

        assert(m_productionResidualEdges[residualPath.from].capacity >= delta);
        m_productionResidualEdges[residualPath.from].capacity -= delta;
        m_productionEdges[residualPath.from].flow += delta;
        assert(m_productionEdges[residualPath.from].flow <= m_productionEdges[residualPath.from].capacity);
        if (m_productionResidualEdges[residualPath.from].capacity == 0) {
            productionSaturated = true;
        }
        for (std::size_t i = node; i < residualPath.from; ++i) {
            assert(m_backwardResidualEdges[i].capacity >= delta);
            m_backwardResidualEdges[i].capacity -= delta;
            m_backwardEdges[i].flow += delta;
            assert(m_backwardEdges[i].flow <= m_backwardEdges[i].capacity);
            m_forwardResidualEdges[i].capacity = m_backwardEdges[i].flow;
            m_forwardResidualEdges[i].cost = -m_backwardEdges[i].cost;
            if (m_backwardResidualEdges[i].capacity == 0) {
                zeroBackwardResidualEdges.push_back(i);
            }
        }

        // Update the backward residual path segments.
        assert(!backwardResidualPathSegments.empty());
        auto& backwardResidualPathSegment = backwardResidualPathSegments.front();

        assert(!backwardResidualPathSegment.first.empty());
        if (productionSaturated || !zeroBackwardResidualEdges.empty()) {
            backwardResidualPathSegment.first.pop_front();
            if (backwardResidualPathSegment.first.empty()) {
                backwardResidualPathSegments.pop_front();
                return;
            }
        }

        if (zeroBackwardResidualEdges.empty()) {
            return;
        }

        // Split the first `backwardResidualPathSegment`.
        std::list<BackwardResidualPathSegment> newBackwardResidualPathSegments;
        std::size_t start = node;
        std::size_t end = start;
        for (std::size_t i: zeroBackwardResidualEdges) {
            end = i;
            if (start < end) {
                newBackwardResidualPathSegments.push_back(std::make_pair(std::list<ResidualPath>(),
                                                                         std::make_pair(start, end)));
            }
            start = i + 1;
            end = start;
        }
        // Last list
        end = backwardResidualPathSegment.second.second;
        if (start < end) {
            newBackwardResidualPathSegments.push_back(std::make_pair(std::list<ResidualPath>(),
                                                                     std::make_pair(start, end)));
        }

        for (const auto& resPath: backwardResidualPathSegment.first) {
            auto it = newBackwardResidualPathSegments.begin();
            while (it != newBackwardResidualPathSegments.end()
                   && (resPath.from < it->second.first || resPath.from > it->second.second)) {
                ++it;
            }

            if (it == newBackwardResidualPathSegments.end()) {
                continue;
            }

            if (resPath.from != it->second.first) {
                it->first.push_back(resPath);
            }
        }

        // Only keep the non-empty ones.
        newBackwardResidualPathSegments.remove_if([](const BackwardResidualPathSegment& residualPathSegment) {
            return residualPathSegment.first.empty();
        });

        // Finally, update the full segments.
        backwardResidualPathSegments.pop_front();
        if (!newBackwardResidualPathSegments.empty()) {
            backwardResidualPathSegments.insert(backwardResidualPathSegments.begin(),
                                                newBackwardResidualPathSegments.begin(),
                                                newBackwardResidualPathSegments.end());
        }
    }
}
