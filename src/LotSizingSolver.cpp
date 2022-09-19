//
//  LotSizingSolver.cpp
//  LotSizingSolver
//
//  Created by ChengLu on 11/24/21.
//  Copyright © 2021 Cheng Lu. All rights reserved.
//

#include "LotSizingSolver.hpp"

#include "dp_array.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <utility>

// MARK: - File static functions

#ifdef DEBUG_LOTSIZING
static void print(const std::list<ResidualPath>& residualPaths, double offset = 0) {
    for (const auto& residualPath: residualPaths) {
        std::cout << "(" << residualPath.from << " of " << residualPath.cost - offset << ") ";
    }
}

static void print(const std::list<BackwardResidualPathSegment>& backwardResidualPathSegments,
                  const std::vector<double>& prefixCostSums) {
    for (const auto& segment: backwardResidualPathSegments) {
        std::cout << "[" << segment.second.first << "," << segment.second.second << "]: ";
        double offset = 0;
        if (segment.second.first > 0) {
            offset = prefixCostSums[segment.second.first - 1];
        }
        print(segment.first, offset);
        std::cout << "\n";
    }
}
#endif

// MARK: - Auxiliary function

const double PW_INC_UNIF_LEFT = 0.1;
const double PW_INC_UNIF_RIGHT = 2.0;

std::random_device rd;
std::mt19937 gen(rd());

std::vector<double> generateConvexCosts(int bkpNum, double leftBkp, double leftSlope) {
    assert(bkpNum > 0);
    // Guarantee bounded minimum.
    assert(leftSlope < 0);
    // Prepare distributions
    std::uniform_real_distribution<double> inc_distribution(PW_INC_UNIF_LEFT, PW_INC_UNIF_RIGHT);
    std::vector<double> pw(bkpNum * 2 + 1, 0);
    pw[0] = leftSlope;
    for (int i = 0; i < bkpNum; ++i) {
        pw[2 * i + 1] = leftBkp;
        double inc = inc_distribution(gen);
        leftSlope += inc;
        // Force the last half pieces to have positive gradients
        // in order to make the problem have bounded minimum.
        if (i >= bkpNum / 2 && leftSlope <= 0) {
            leftSlope = 1;
        }
        pw[2 * i + 2] = leftSlope;
        inc = inc_distribution(gen);
        leftBkp += inc;
    }

    return pw;
}

// MARK: - Convex/Linear edge generation

ConvexEdge::ConvexEdge(const std::vector<double>& costs) {
    assert(costs.size() % 2 == 1);
    for (int i = 0; i < costs.size() / 2; ++i) {
        _slopes.push_back(costs[i * 2]);
        _bkps.push_back(costs[i * 2 + 1]);
    }
    _slopes.push_back(costs[costs.size() - 1]);
    _bkpValues.push_back(0);
    for (int i = 1; i < _bkps.size(); ++i) {
        double segValue = (_bkps[i] - _bkps[i - 1]) * _slopes[i];
        _bkpValues.push_back(_bkpValues[i - 1] + segValue);
    }
}

double ConvexEdge::cost(double flow, double delta) {
    assert(delta > 0);
    return (_cost(flow + delta) - _cost(flow)) / delta;
}

double ConvexEdge::_cost(double flow) {
    // Binary search to find the segment.
    std::size_t head = 0;
    std::size_t tail = _bkps.size();
    while (head < tail) {
        auto index = (head + tail) / 2;
        if (flow <= _bkps[index]) {
            if (index == head || flow > _bkps[index - 1]) {
                return _bkpValues[index] + (flow - _bkps[index]) * _slopes[index];
            }
            tail = index;
        } else {
            if (index + 1 == tail || flow <= _bkps[index + 1]) {
                return _bkpValues[index] + (flow - _bkps[index]) * _slopes[index + 1];
            }
            head = index + 1;
        }
    }

    // Expecte
    assert(false);
    return -1;
}

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
    m_productionResidualEdges.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        m_demands[i] = demands[i];
        totalDemand += demands[i];

        auto edge = std::make_shared<Edge>(productionCosts[i]);
        m_productionEdges[i] = edge;
        edge->capacity = productionCapacities[i];
        edge->flow = 0;
        assert(edge->cost() >= 0);
        totalProductionCapacity += edge->capacity;
        // Partial-sum satisfactory as well.
        assert(totalDemand <= totalProductionCapacity);
        m_costSum += edge->cost();

        // Production residual edges
        edge = std::make_shared<Edge>(productionCosts[i]);
        m_productionResidualEdges[i] = edge;
        edge->capacity = productionCapacities[i];
        edge->flow = 0;
    }

    // Inventory forward edges
    m_forwardEdges.resize(n - 1);
    m_forwardResidualEdges.resize(n - 1);
    for (std::size_t i = 0; i < n - 1; ++i) {
        auto edge = std::make_shared<Edge>(forwardCosts[i]);
        m_forwardEdges[i] = edge;
        edge->capacity = forwardCapacities[i];
        edge->flow = 0;
        assert(edge->cost() >= 0);
        m_costSum += edge->cost();

        // Forward residual edges
        edge = std::make_shared<Edge>(forwardCosts[i]);
        m_forwardResidualEdges[i] = edge;
        edge->capacity = forwardCapacities[i];
        edge->flow = 0;
    }

    m_costBound = m_costSum * 10;
}

double ForwardGraph::cost() const {
    double cost = 0;
    for (const auto& edge: m_productionEdges) {
        cost += edge->cost() * edge->flow;
    }

    for (const auto& edge: m_forwardEdges) {
        cost += edge->cost() * edge->flow;
    }

    return cost;
}

bool ForwardGraph::solve() {
    std::list<ResidualPath> residualPaths;

    for (std::size_t i = 0; i < m_n; ++i) {
        elongateAndAdd(i, residualPaths);
        uint32_t demand = m_demands[i];
        while (demand > 0) {
            if (!augmentAndUpdate(i, residualPaths, demand)) {
                return false;
            }
        }
    }

    return true;
}

bool ForwardGraph::fasterSolve() {
    // Initiate the capacity dynamic path and the residual path cost dynamic path.
    std::vector<uint32_t> capacities;
    for (std::size_t i = 0; i < m_n - 1; ++i) {
        capacities.push_back(m_forwardResidualEdges[i]->capacity);
    }
    dp_array<uint32_t> capacityDP(capacities);

    std::vector<double> costs(m_n, 0);
    for (std::size_t i = 0; i < m_n; ++i) {
        if (m_productionResidualEdges[i]->capacity > 0) {
            costs[i] = m_productionResidualEdges[i]->cost();
        } else {
            costs[i] = m_costBound;
        }
    }
    dp_array<double> costDP(costs);

    dp_array<double> flowDP(std::vector<double>(m_n - 1, 0));

    std::size_t start = 0;
    for (std::size_t i = 0; i < m_n; ++i) {
        fastElongateAndAdd(i, start, capacityDP, costDP);
        uint32_t demand = m_demands[i];
        while (demand > 0) {
            if (!fastAugmentAndUpdate(i, start,
                                      capacityDP, costDP, flowDP,
                                      demand)) {
                return false;
            }
        }
    }
    std::vector<double> flows;
    assert(flowDP.vectorize(flows));
    for (std::size_t i = 0; i < m_n - 1; ++i) {
        m_forwardEdges[i]->flow = flows[i];
    }

    return true;
}

bool ForwardGraph::constraintsSatisfied() const {
    // Capacity constraints
    for (const auto& edge: m_productionEdges) {
        if (edge->flow > edge->capacity) {
            return false;
        }
    }

    for (const auto& edge: m_forwardEdges) {
        if (edge->flow > edge->capacity) {
            return false;
        }
    }

    // Flow balance constraints
    // First node
    if (static_cast<uint32_t>(m_productionEdges[0]->flow - m_forwardEdges[0]->flow) != m_demands[0]) {
        return false;
    }

    for (std::size_t i = 1; i < m_n - 1; ++i) {
        if (static_cast<uint32_t>(m_productionEdges[i]->flow + m_forwardEdges[i - 1]->flow - m_forwardEdges[i]->flow) != m_demands[i]) {
            return false;
        }
    }

    if (static_cast<uint32_t>(m_productionEdges[m_n - 1]->flow + m_forwardEdges[m_n - 2]->flow) != m_demands[m_n - 1]) {
        return false;
    }

    return true;
}

bool ForwardGraph::isOptimal() const {
    std::vector<double> d(m_n, std::numeric_limits<double>::max());

    // Bellman-Ford shortest path algorithm.

    for (std::size_t i = 0; i < m_n; ++i) {
        // Production arcs
        for (std::size_t j = 0; j < m_n; ++j) {
            const auto& edge = m_productionEdges[j];
            if (edge->flow < edge->capacity) {
                if (edge->cost() < d[j]) {
                    d[j] = edge->cost();
                }
            }

            if (edge->flow > 0) {
                if (d[j] - edge->cost() < 0) {
                    return false;
                }
            }
        }

        for (std::size_t j = 0; j < m_n - 1; ++j) {
            const auto& edge = m_forwardEdges[j];
            if (edge->flow < edge->capacity) {
                if (d[j] + edge->cost() < d[j + 1]) {
                    d[j + 1] = d[j] + edge->cost();
                }
            }

            if (edge->flow > 0) {
                if (d[j + 1] - edge->cost() < d[j]) {
                    d[j] = d[j + 1] - edge->cost();
                }
            }
        }
    }

    // Verify no-negative cycles.
    for (std::size_t i = 0; i < m_n; ++i) {
        const auto& edge = m_productionEdges[i];
        if (edge->flow < edge->capacity) {
            if (edge->cost() < d[i]) {
                return false;
            }
        }

        if (edge->flow > 0) {
            if (d[i] - edge->cost() < 0) {
                return false;
            }
        }
    }

    for (std::size_t i = 0; i < m_n - 1; ++i) {
        const auto& edge = m_forwardEdges[i];
        if (edge->flow < edge->capacity) {
            if (d[i] + edge->cost() < d[i + 1]) {
                return false;
            }
        }

        if (edge->flow > 0) {
            if (d[i + 1] - edge->cost() < d[i]) {
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

bool ForwardGraph::isInfinite(double value) {
    return value > m_costSum;
}

#ifdef DEBUG_LOTSIZING
void ForwardGraph::print() const {
    // Production edges.
    for (const auto& edge: m_productionEdges) {
        std::cout << "(" << edge->capacity << "," << edge->cost() << "," << edge->flow << ")";
    }
    std::cout << "\n";

    // Inventory forward edges.
    for (const auto& edge: m_forwardEdges) {
        std::cout << "(" << edge->capacity << "," << edge->cost() << "," << edge->flow << ")";
    }
    std::cout << "\n";
}

void ForwardGraph::printResiduals() const {
    // Production residuals
    for (const auto& edge: m_productionResidualEdges) {
        std::cout << "(" << edge->capacity << "," << edge->cost() << ")";
    }
    std::cout << "\n";

    // Forward residuals
    for (const auto& edge: m_forwardResidualEdges) {
        std::cout << "(" << edge->capacity << "," << edge->cost() << ")";
    }
    std::cout << "\n";
}
#endif

// MARK: - Private

void ForwardGraph::elongateAndAdd(std::size_t node, std::list<ResidualPath>& residualPaths) {
    if (node > 0) {
        // First append arc (node - 1, node) to each existing residual path.
        const auto& residualEdge = m_forwardResidualEdges[node - 1];

        if (residualEdge->capacity > 0) {
            for (auto it = residualPaths.begin(); it != residualPaths.end(); ++it) {
                it->cost += residualEdge->cost();
            }
        } else {
            // Clear all existing residual paths.
            residualPaths.clear();
        }
    }

    // Insert the residual path (0) -> (node)
    if (m_productionResidualEdges[node]->capacity > 0) {
        ResidualPath residualPath;
        residualPath.from = node;
        residualPath.cost = m_productionResidualEdges[node]->cost();
        // Insert
        orderedInsert(residualPath, residualPaths);
    }
}

bool ForwardGraph::augmentAndUpdate(std::size_t node, std::list<ResidualPath>& residualPaths, uint32_t& demand) {
    assert(demand > 0);

    if (residualPaths.empty()) {
        return false;
    }

    // Get the minimum cost residual path.
    const auto& residualPath = residualPaths.front();

    // Find residual capacity
    uint32_t residualCapacity = m_productionResidualEdges[residualPath.from]->capacity;
    for (std::size_t i = residualPath.from; i < node; ++i) {
        residualCapacity = std::min(residualCapacity, m_forwardResidualEdges[i]->capacity);
    }
    assert(residualCapacity > 0);

    auto delta = std::min(demand, residualCapacity);
    assert(delta > 0);

    // Update demand
    demand -= delta;

    // Update flows and residual edges.
    bool productionSaturated = false;
    std::size_t lastZeroResCapNode = m_n;
    assert(m_productionResidualEdges[residualPath.from]->capacity >= delta);
    m_productionResidualEdges[residualPath.from]->capacity -= delta;
    m_productionEdges[residualPath.from]->flow += delta;
    if (m_productionResidualEdges[residualPath.from]->capacity == 0) {
        productionSaturated = true;
    }
    for (std::size_t i = residualPath.from; i < node; ++i) {
        assert(m_forwardResidualEdges[i]->capacity >= delta);
        m_forwardResidualEdges[i]->capacity -= delta;
        // For forward graph, the forward edges are the same as the forward residual edges.
        m_forwardEdges[i]->flow += delta;
        if (m_forwardResidualEdges[i]->capacity == 0) {
            lastZeroResCapNode = i;
        }
    }

    // Update the residual paths.
    if (productionSaturated) {
        residualPaths.pop_front();
    }

    if (residualPath.from == node) {
        return true;
    }

    if (lastZeroResCapNode < m_n) {
        residualPaths.remove_if([&node, &lastZeroResCapNode](const ResidualPath& residualPath) {
            return residualPath.from != node && residualPath.from <= lastZeroResCapNode;
        });
    }

    return true;
}

void ForwardGraph::fastElongateAndAdd(std::size_t node, std::size_t& start,
                                      dp_array<uint32_t>& capacityDP, dp_array<double>& costDP) {
    if (node > 0) {
        auto nodeInt = static_cast<int>(node);
        auto startInt = static_cast<int>(start);
        auto capacity = capacityDP.edge_cost(nodeInt - 1);
        const auto& residualEdge = m_forwardResidualEdges[node - 1];
        if (capacity > 0) {
            costDP.update_constant(startInt, nodeInt, residualEdge->cost());
        } else {
            start = node;
        }
    }
}

bool ForwardGraph::fastAugmentAndUpdate(std::size_t node, std::size_t& start,
                                        dp_array<uint32_t>& capacityDP, dp_array<double>& costDP, dp_array<double>& flowDP,
                                        uint32_t& demand) {
    assert(demand > 0);

    auto nodeInt = static_cast<int>(node);
    auto startInt = static_cast<int>(start);

    // Get the minimum cost residual path.
    int from = 0;
    auto minCost = costDP.min_cost_last(startInt, nodeInt + 1, from);
    assert(from >= 0);

    if (!minCost.has_value() || isInfinite(*minCost)) {
        return false;
    }

    int tmp;
    uint32_t capacity = m_productionResidualEdges[from]->capacity;
    if (from < node) {
        capacity = std::min(capacity, *(capacityDP.min_cost_last(from, nodeInt, tmp)));
    }
    assert(capacity > 0);

    auto delta = std::min(demand, capacity);
    assert(delta > 0);

    // Update demand
    demand -= delta;

    // Update flows and residual edges.
    bool productionSaturated = false;
    assert(m_productionResidualEdges[from]->capacity >= delta);
    m_productionResidualEdges[from]->capacity -= delta;
    m_productionEdges[from]->flow += delta;
    if (m_productionResidualEdges[from]->capacity == 0) {
        productionSaturated = true;
    }

    capacityDP.update_constant(from, nodeInt, -delta);
    flowDP.update_constant(from, nodeInt, delta);

    if (productionSaturated) {
        // Essentially remove this residual path.
        costDP.update_constant(from, from + 1, m_costBound);
    }

    if (from == node) {
        return true;
    }

    capacity = *capacityDP.min_cost_last(from, nodeInt, tmp);
    if (capacity == 0) {
        start = tmp + 1;
    }

    return true;
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
    m_backwardResidualEdges.resize(n - 1);
    m_backwardResidualEdgeCostSums.resize(n - 1);
    for (std::size_t i = 0; i < n - 1; ++i) {
        auto edge = std::make_shared<Edge>(backwardCosts[i]);
        m_backwardEdges[i] = edge;
        edge->capacity = backwardCapacities[i];
        edge->flow = 0;
        m_costSum += edge->cost();
        assert(edge->cost() >= 0);
        if (i == 0) {
            m_backwardResidualEdgeCostSums[i] = edge->cost();
        } else {
            m_backwardResidualEdgeCostSums[i] = m_backwardResidualEdgeCostSums[i - 1] + edge->cost();
        }

        // Backward residual edges
        edge = std::make_shared<Edge>(backwardCosts[i]);
        m_backwardResidualEdges[i] = edge;
        edge->capacity = backwardCapacities[i];
        edge->flow = 0;
    }

    m_costBound = m_costSum * 10;

#ifdef DEBUG_LOTSIZING
    // Print demands
    std::cout << "Start a new forward-backward graph with demands:\n";
    for (uint32_t demand: m_demands) {
        std::cout << demand << ",";
    }
    std::cout << "\n";
#endif
}

double ForwardBackwardGraph::cost() const {
    double cost = ForwardGraph::cost();

    for (const auto& edge: m_backwardEdges) {
        cost += edge->cost() * edge->flow;
    }

    return cost;
}

bool ForwardBackwardGraph::solve() {
    std::list<ResidualPath> forwardResidualPaths;
    std::list<BackwardResidualPathSegment> backwardResidualPathSegments;

    // Initialize the forward and backward residual paths.
    if (m_productionResidualEdges[0]->capacity > 0) {
        ResidualPath residualPath;
        residualPath.from = 0;
        residualPath.cost = m_productionResidualEdges[0]->cost();
        orderedInsert(residualPath, forwardResidualPaths);
    }

    // Accumulated backward residual cost.
    std::list<ResidualPath> backwardResidualPath;
    std::size_t start = 0;
    std::size_t end = start;
    for (std::size_t i = 1; i < m_n; ++i) {
        if (m_backwardResidualEdges[i - 1]->capacity == 0) {
            if (start < end && !backwardResidualPath.empty()) {
                backwardResidualPathSegments.push_back(std::make_pair(backwardResidualPath,
                                                                      std::make_pair(start, end)));
            }
            start = i;
            end = start;
            backwardResidualPath.clear();
        } else {
            // Ordered insert a backward residual path.
            if (m_productionResidualEdges[i]->capacity > 0 && i != start) {
                ResidualPath residualPath;
                residualPath.from = i;
                residualPath.cost = m_backwardResidualEdgeCostSums[i - 1] + m_productionResidualEdges[i]->cost();
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

#ifdef DEBUG_LOTSIZING
    // Initial status
    std::cout << "Initial:\n";
    print(forwardResidualPaths, backwardResidualPathSegments);
#endif

    for (std::size_t i = 0; i < m_n; ++i) {
        uint32_t demand = m_demands[i];
#ifdef DEBUG_LOTSIZING
        std::cout << "Start to push for Node " << i << ":\n";
#endif
        while (demand > 0) {
            if (!augmentAndUpdate(i, forwardResidualPaths, backwardResidualPathSegments, demand)) {
                return false;
            }
#ifdef DEBUG_LOTSIZING
            print(forwardResidualPaths, backwardResidualPathSegments);
            std::cout << "Push down demand to " << demand << std::endl;
#endif
        }
#ifdef DEBUG_LOTSIZING
        std::cout << "Done pushing for Node " << i << ".\n";
#endif
        if (i < m_n - 1) {
            elongateAndUpdate(i + 1, forwardResidualPaths, backwardResidualPathSegments);
#ifdef DEBUG_LOTSIZING
            print(forwardResidualPaths, backwardResidualPathSegments);
            std::cout << "Done for preparing for Node " << i + 1 << ".\n";
#endif
        }
    }

#ifdef DEBUG_LOTSIZING
    // Illustrate the residual edges.
    std::cout << "Final status >>>>>>>>>>\n";
    print(forwardResidualPaths, backwardResidualPathSegments);
    std::cout << "Residual edges:\n";
    printResiduals();
#endif

    return true;
}

bool ForwardBackwardGraph::fasterSolve() {
    // Initiate the dynamic paths.
    std::vector<uint32_t> forwardResCapacities(m_n - 1,0);
    dp_array forwardResCapacityDP(forwardResCapacities);

    std::vector<uint32_t> backwardResCapacities;
    for (std::size_t i = 0; i < m_n - 1; ++i) {
        backwardResCapacities.push_back(m_backwardResidualEdges[i]->capacity);
    }
    dp_array backwardResCapacityDP(backwardResCapacities);

    std::vector<double> costs(m_n, 0);
    // Node 0 -- Forward
    if (m_productionResidualEdges[0]->capacity > 0) {
        costs[0] = m_productionResidualEdges[0]->cost();
    } else {
        costs[0] = m_costBound;
    }

    // Remaining nodes -- Backward
    for (std::size_t i = 1; i < m_n; ++i) {
        if (m_productionResidualEdges[i]->capacity > 0) {
            costs[i] = m_productionResidualEdges[i]->cost() + m_backwardResidualEdgeCostSums[i - 1];
        } else {
            costs[i] = m_costBound;
        }
    }
    dp_array costDP(costs);

    dp_array forwardResFlowDP(std::vector<double>(m_n - 1, 0));
    dp_array backwardResFlowDP(std::vector<double>(m_n - 1, 0));

    // Find backward residual segments.
    std::list<SegmentRange> backwardResidualSegments;
    std::size_t start = 0;
    std::size_t end = start;
    for (std::size_t i = 1; i < m_n; ++i) {
        if (m_backwardResidualEdges[i - 1]->capacity == 0) {
            if (start < end) {
                backwardResidualSegments.push_back(SegmentRange(start, end));
            }
            start = i;
            end = start;
        } else {
            end = i;
        }
    }

    // Last segment
    assert(end == m_n - 1);
    if (start < end) {
        backwardResidualSegments.push_back(SegmentRange(start, end));
    }

    start = 0;
    for (std::size_t i = 0; i < m_n; ++i) {
        uint32_t demand = m_demands[i];
        while (demand > 0) {
            if (!fastAugmentAndUpdate(i, start, backwardResidualSegments,
                                      forwardResCapacityDP, backwardResCapacityDP,
                                      costDP, forwardResFlowDP, backwardResFlowDP,
                                      demand)) {
                return false;
            }
        }

        if (i < m_n - 1) {
            fastElongateAndUpdate(i + 1, forwardResFlowDP, backwardResFlowDP,
                                  start, backwardResidualSegments,
                                  forwardResCapacityDP, backwardResCapacityDP, costDP);
        }
    }

    mergeFlowSolutions(forwardResFlowDP, backwardResFlowDP);

    return true;
}

bool ForwardBackwardGraph::constraintsSatisfied() const {
    // Capacity constraints
    for (const auto& edge: m_productionEdges) {
        if (edge->flow > edge->capacity) {
            return false;
        }
    }

    for (const auto& edge: m_forwardEdges) {
        if (edge->flow > edge->capacity) {
            return false;
        }
    }

    for (const auto& edge: m_backwardEdges) {
        if (edge->flow > edge->capacity) {
            return false;
        }
    }

    for (std::size_t i = 0; i < m_n - 1; ++i) {
        if (m_forwardEdges[i]->flow > 0 && m_backwardEdges[i]->flow > 0) {
            return false;
        }
    }

    // Flow balance constraints
    // First node
    if (static_cast<uint32_t>(m_productionEdges[0]->flow + m_backwardEdges[0]->flow - m_forwardEdges[0]->flow)
        != m_demands[0]) {
        return false;
    }

    for (std::size_t i = 1; i < m_n - 1; ++i) {
        if (static_cast<uint32_t>(m_productionEdges[i]->flow + m_forwardEdges[i - 1]->flow + m_backwardEdges[i]->flow
            - m_forwardEdges[i]->flow - m_backwardEdges[i - 1]->flow)
            != m_demands[i]) {
            return false;
        }
    }

    if (static_cast<uint32_t>(m_productionEdges[m_n - 1]->flow + m_forwardEdges[m_n - 2]->flow - m_backwardEdges[m_n - 2]->flow)
        != m_demands[m_n - 1]) {
        return false;
    }

    return true;
}

const double kEps = 1e-6;

static double truncatedValue(double value) {
    if (std::fabs(value) < kEps) {
        return 0;
    }
    return value;
}

bool ForwardBackwardGraph::isOptimal() const {
    std::vector<double> d(m_n, std::numeric_limits<double>::max());

    // Bellman-Ford shortest path algorithm.

    for (std::size_t i = 0; i < m_n; ++i) {
        // Production arcs
        for (std::size_t j = 0; j < m_n; ++j) {
            const auto& edge = m_productionEdges[j];
            if (edge->flow < edge->capacity) {
                double value = truncatedValue(edge->cost());
                if (value < d[j]) {
                    d[j] = value;
                }
            }

            if (edge->flow > 0) {
                double value = truncatedValue(d[j] - edge->cost());
                if (value < 0) {
                    return false;
                }
            }
        }

        for (std::size_t j = 0; j < m_n - 1; ++j) {
            const auto& fEdge = m_forwardEdges[j];
            if (fEdge->flow < fEdge->capacity) {
                double value = truncatedValue(d[j] + fEdge->cost());
                if (value < d[j + 1]) {
                    d[j + 1] = value;
                }
            }
            if (fEdge->flow > 0) {
                double value = truncatedValue(d[j + 1] - fEdge->cost());
                if (value < d[j]) {
                    d[j] = value;
                }
            }

            const auto& bEdge = m_backwardEdges[j];
            if (bEdge->flow < bEdge->capacity) {
                double value = truncatedValue(d[j + 1] + bEdge->cost());
                if (value < d[j]) {
                    d[j] = value;
                }
            }

            if (bEdge->flow > 0) {
                double value = d[j] - bEdge->cost();
                if (value < d[j + 1]) {
                    d[j + 1] = value;
                }
            }
        }
    }

    // Verify no-negative cycles.
    for (std::size_t i = 0; i < m_n; ++i) {
        const auto& edge = m_productionEdges[i];
        if (edge->flow < edge->capacity) {
            if (truncatedValue(edge->cost()) + kEps < d[i]) {
                return false;
            }
        }

        if (edge->flow > 0) {
            if (truncatedValue(d[i] - edge->cost()) + kEps < 0) {
                return false;
            }
        }
    }

    for (std::size_t i = 0; i < m_n - 1; ++i) {
        const auto& fEdge = m_forwardEdges[i];
        if (fEdge->flow < fEdge->capacity) {
            if (truncatedValue(d[i] + fEdge->cost()) + kEps < d[i + 1]) {
                return false;
            }
        }

        if (fEdge->flow > 0) {
            if (truncatedValue(d[i + 1] - fEdge->cost()) + kEps < d[i]) {
                return false;
            }
        }

        const auto& bEdge = m_backwardEdges[i];
        if (bEdge->flow < bEdge->capacity) {
            if (truncatedValue(d[i + 1] + bEdge->cost()) + kEps < d[i]) {
                return false;
            }
        }

        if (bEdge->flow > 0) {
            if (truncatedValue(d[i] - bEdge->cost()) + kEps < d[i + 1]) {
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

    if (residualEdge->capacity > 0) {
        for (auto it = forwardResidualPaths.begin(); it != forwardResidualPaths.end(); ++it) {
            it->cost += residualEdge->cost();
        }
    } else {
        // Clear all existing forward residual paths.
        forwardResidualPaths.clear();
    }

    // Insert the residual path (0) -> (node)
    if (m_productionResidualEdges[node]->capacity > 0) {
        ResidualPath residualPath;
        residualPath.from = node;
        residualPath.cost = m_productionResidualEdges[node]->cost();
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

bool ForwardBackwardGraph::augmentAndUpdate(std::size_t node, std::list<ResidualPath>& forwardResidualPaths,
                                            std::list<BackwardResidualPathSegment>& backwardResidualPathSegments,
                                            uint32_t& demand) {
    assert(demand > 0);

    if (forwardResidualPaths.empty() && backwardResidualPathSegments.empty()) {
        return false;
    }

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
    uint32_t residualCapacity = m_productionResidualEdges[residualPath.from]->capacity;

    if (node >= residualPath.from) {
        // Forward residual path
        for (std::size_t i = residualPath.from; i < node; ++i) {
            residualCapacity = std::min(residualCapacity, m_forwardResidualEdges[i]->capacity);
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

        assert(m_productionResidualEdges[residualPath.from]->capacity >= delta);
        m_productionResidualEdges[residualPath.from]->capacity -= delta;
        m_productionEdges[residualPath.from]->flow += delta;
        assert(m_productionEdges[residualPath.from]->flow <= m_productionEdges[residualPath.from]->capacity);
        if (m_productionResidualEdges[residualPath.from]->capacity == 0) {
            producationSaturated = true;
        }
        for (std::size_t i = residualPath.from; i < node; ++i) {
            assert(m_forwardResidualEdges[i]->capacity >= delta);
            m_forwardResidualEdges[i]->capacity -= delta;
            // Two cases
            if (m_backwardEdges[i]->flow > 0) {
                // Reverse arc
                assert(m_backwardEdges[i]->flow >= delta);
                assert(static_cast<uint32_t>(m_forwardEdges[i]->flow) == 0);
                m_backwardEdges[i]->flow -= delta;

                m_backwardResidualEdges[i]->capacity = m_backwardEdges[i]->capacity - m_backwardEdges[i]->flow;
                if (m_forwardResidualEdges[i]->capacity == 0) {
                    assert(static_cast<uint32_t>(m_backwardEdges[i]->flow) == 0);
                    // Turn to the regular arc
                    if (m_forwardEdges[i]->capacity > 0) {
                        double costDelta = m_backwardEdges[i]->cost() + m_forwardEdges[i]->cost();
                        residualCostUpdatedEdges.push_back(std::make_pair(i, costDelta));
                        m_forwardResidualEdges[i] = std::make_shared<Edge>(m_forwardEdges[i]->cost());
                        m_forwardResidualEdges[i]->capacity = m_forwardEdges[i]->capacity;
                    } else {
                        lastZeroResCapNode = i;
                    }
                }
            } else {
                // Regular arc
                m_forwardEdges[i]->flow += delta;
                m_backwardResidualEdges[i] = std::make_shared<Edge>(-m_forwardEdges[i]->cost());
                m_backwardResidualEdges[i]->capacity = m_forwardEdges[i]->flow;
                assert(m_forwardEdges[i]->flow <= m_forwardEdges[i]->capacity);
                if (m_forwardResidualEdges[i]->capacity == 0) {
                    lastZeroResCapNode = i;
                }
            }
        }

        // Update the forward residual paths.
        if (producationSaturated || lastZeroResCapNode < m_n) {
            forwardResidualPaths.pop_front();
        }

        if (residualPath.from == node) {
            return true;
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
            residualCapacity = std::min(residualCapacity, m_backwardResidualEdges[i]->capacity);
        }
        assert(residualCapacity > 0);

        auto delta = std::min(demand, residualCapacity);
        assert(delta > 0);

        // Update demand
        demand -= delta;

        // Update flows and residual edges
        std::vector<std::size_t> zeroBackwardResidualEdges;
        bool productionSaturated = false;

        assert(m_productionResidualEdges[residualPath.from]->capacity >= delta);
        m_productionResidualEdges[residualPath.from]->capacity -= delta;
        m_productionEdges[residualPath.from]->flow += delta;
        assert(m_productionEdges[residualPath.from]->flow <= m_productionEdges[residualPath.from]->capacity);
        if (m_productionResidualEdges[residualPath.from]->capacity == 0) {
            productionSaturated = true;
        }
        for (std::size_t i = node; i < residualPath.from; ++i) {
            assert(m_backwardResidualEdges[i]->capacity >= delta);
            m_backwardResidualEdges[i]->capacity -= delta;
            m_backwardEdges[i]->flow += delta;
            assert(m_backwardEdges[i]->flow <= m_backwardEdges[i]->capacity);
            m_forwardResidualEdges[i] = std::make_shared<Edge>(-m_backwardEdges[i]->cost());
            m_forwardResidualEdges[i]->capacity = m_backwardEdges[i]->flow;
            if (m_backwardResidualEdges[i]->capacity == 0) {
                zeroBackwardResidualEdges.push_back(i);
            }
        }

        // Update the backward residual path segments.
        assert(!backwardResidualPathSegments.empty());
        auto& backwardResidualPathSegment = backwardResidualPathSegments.front();

        assert(!backwardResidualPathSegment.first.empty());
        if (productionSaturated) {
            backwardResidualPathSegment.first.pop_front();
            if (backwardResidualPathSegment.first.empty()) {
                backwardResidualPathSegments.pop_front();
                return true;
            }
        }

        if (zeroBackwardResidualEdges.empty()) {
            return true;
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

    return true;
}

void ForwardBackwardGraph::fastElongateAndUpdate(std::size_t node, const dp_array<double>& forwardResFlowDP, const dp_array<double>& backwardResFlowDP,
                                                 std::size_t& start, std::list<SegmentRange>& backwardResidualSegments,
                                                 dp_array<uint32_t>& forwardResCapacityDP, dp_array<uint32_t>& backwardResCapacityDP,
                                                 dp_array<double>& costDP) {
    assert(node > 0);
    auto nodeInt = static_cast<int>(node);
    auto startInt = static_cast<int>(start);

    // MARK: - Forward residual paths.

    // Note: Check backward residual flow first.

    auto backwardResFlow = *backwardResFlowDP.edge_cost(nodeInt - 1);

    if (backwardResFlow > 0) {
        // Use backward flow as forward residual capacity.
        auto capacity = backwardResFlow;
        double cost = -m_backwardResidualEdges[nodeInt - 1]->cost();
        forwardResCapacityDP.update_constant(nodeInt - 1, nodeInt, capacity);
        costDP.update_constant(startInt, nodeInt, cost);
    } else {
        // Use original forward residual capacity.
        auto capacity = m_forwardResidualEdges[node - 1]->capacity;
        auto cost = m_forwardResidualEdges[node - 1]->cost();
        if (capacity > 0) {
            forwardResCapacityDP.update_constant(nodeInt - 1, nodeInt, capacity);
            costDP.update_constant(startInt, nodeInt, cost);
        } else {
            start = node;
        }
    }

    costDP.update_constant(nodeInt, nodeInt + 1, -m_backwardResidualEdgeCostSums[node - 1]);

    // MARK: - Backward residual paths.

    if (backwardResidualSegments.empty()) {
        return;
    }

    auto& backwardResidualSegment = backwardResidualSegments.front();
    if (backwardResidualSegment.first > node) {
        return;
    }

    assert(backwardResidualSegment.first >= node - 1
           && backwardResidualSegment.first <= node);
    backwardResidualSegment.first = node;
    if (backwardResidualSegment.first >= backwardResidualSegment.second) {
        backwardResidualSegments.pop_front();
        return;
    }
}

bool ForwardBackwardGraph::fastAugmentAndUpdate(std::size_t node, std::size_t& start, std::list<SegmentRange>& backwardResidualSegments,
                                                dp_array<uint32_t>& forwardResCapacityDP, dp_array<uint32_t>& backwardResCapacityDP,
                                                dp_array<double>& costDP,
                                                dp_array<double>& forwardResFlowDP, dp_array<double>& backwardResFlowDP,
                                                uint32_t& demand) {
    assert(demand > 0);

    auto nodeInt = static_cast<int>(node);
    auto startInt = static_cast<int>(start);

    // Get the minimum cost out of forward and backward residual paths.
    int from = 0;
    auto minCost = *costDP.min_cost_last(startInt, nodeInt + 1, from);
    assert(from >= 0);
    bool forwardFeasible = !isInfinite(minCost);

    // No feasible forward or backward residual paths.
    if (!forwardFeasible && backwardResidualSegments.empty()) {
        return false;
    }

    if (!backwardResidualSegments.empty()) {
        const auto& backwardResidualSegment = backwardResidualSegments.front();
        const auto firstNode = backwardResidualSegment.first;
        const auto secondNode = backwardResidualSegment.second;
        assert(node <= firstNode);
        if (node == firstNode) {
            int bFrom = 0;
            auto bMinCost = *costDP.min_cost_first(nodeInt + 1, static_cast<int>(secondNode + 1), bFrom);
            if (firstNode > 0) {
                bMinCost -= m_backwardResidualEdgeCostSums[firstNode - 1];
            }

            if (isInfinite(bMinCost) && !forwardFeasible) {
                return false;
            } else if (bMinCost < minCost) {
                from = bFrom;
            }
        }
    }

    // Find residual capacity
    uint32_t residualCapacity = m_productionResidualEdges[from]->capacity;

    int tmp;
    if (nodeInt == from) {
        // Production arc only, simple.
        auto delta = std::min(demand, residualCapacity);
        assert(delta > 0);

        // Update demand
        demand -= delta;

        assert(m_productionResidualEdges[from]->capacity >= delta);
        m_productionResidualEdges[from]->capacity -= delta;
        m_productionEdges[from]->flow += delta;
        if (m_productionResidualEdges[from]->capacity == 0) {
            costDP.update_constant(from, from + 1, m_costBound);
        }
    } else if (nodeInt > from) {
        // Forward residual path
        auto minForwardResCapacity = forwardResCapacityDP.min_cost_last(from, nodeInt, tmp);
        assert(minForwardResCapacity.has_value());
        residualCapacity = std::min(residualCapacity, *minForwardResCapacity);

        assert(residualCapacity > 0);

        auto delta = std::min(demand, residualCapacity);
        assert(delta > 0);

        // Update demand
        demand -= delta;

        // Update flows and residual edges
        assert(m_productionResidualEdges[from]->capacity >= delta);
        m_productionResidualEdges[from]->capacity -= delta;
        m_productionEdges[from]->flow += delta;
        if (m_productionResidualEdges[from]->capacity == 0) {
            costDP.update_constant(from, from + 1, m_costBound);
        }

        forwardResCapacityDP.update_constant(from, nodeInt, -delta);
        forwardResFlowDP.update_constant(from, nodeInt, delta);

        // Check for 0-capacity forward residual arc
        minForwardResCapacity = forwardResCapacityDP.min_cost_last(from, nodeInt, tmp);
        while (minForwardResCapacity.has_value() && (*minForwardResCapacity) == 0) {
            auto forwardResFlow = *forwardResFlowDP.edge_cost(tmp);
            auto backwardResFlow = *backwardResFlowDP.edge_cost(tmp);
            if (backwardResFlow >= forwardResFlow) {
                assert(forwardResFlow == backwardResFlow);
                if (m_forwardResidualEdges[tmp]->capacity > 0) {
                    // Change backward to forward arc.
                    auto costDelta = m_forwardResidualEdges[tmp]->cost() + m_backwardResidualEdges[tmp]->cost();
                    costDP.update_constant(startInt, tmp + 1, costDelta);
                    forwardResCapacityDP.update_constant(tmp, tmp + 1, m_forwardResidualEdges[tmp]->capacity);
                } else {
                    start = static_cast<std::size_t>(tmp + 1);
                    break;
                }
            } else {
                assert(forwardResFlow - backwardResFlow == m_forwardEdges[tmp]->capacity);
                start = static_cast<std::size_t>(tmp + 1);
                break;
            }
            minForwardResCapacity = forwardResCapacityDP.min_cost_last(from, tmp, tmp);
        }
    } else {
        // Backward residual path
        assert(!backwardResidualSegments.empty());

        auto minBackwardResCapacity = backwardResCapacityDP.min_cost_first(nodeInt, from, tmp);
        assert(minBackwardResCapacity.has_value());
        residualCapacity = std::min(residualCapacity, *minBackwardResCapacity);

        assert(residualCapacity > 0);

        auto delta = std::min(demand, residualCapacity);
        assert(delta > 0);

        // Update demand
        demand -= delta;

        // Update flows and residual edges
        assert(m_productionResidualEdges[from]->capacity >= delta);
        m_productionResidualEdges[from]->capacity -= delta;
        m_productionEdges[from]->flow += delta;
        if (m_productionResidualEdges[from]->capacity == 0) {
            costDP.update_constant(from, from + 1, m_costBound);
        }

        backwardResCapacityDP.update_constant(nodeInt, from, -delta);
        backwardResFlowDP.update_constant(nodeInt, from, delta);

        // Check for 0-capacity backward residual arc
        std::list<SegmentRange> splitSegments;
        bool splitHappened = false;
        minBackwardResCapacity = backwardResCapacityDP.min_cost_first(nodeInt, from, tmp);
        int segmentStart = nodeInt;
        while (minBackwardResCapacity.has_value() && (*minBackwardResCapacity) == 0) {
            splitHappened = true;
            if (segmentStart < tmp) {
                splitSegments.push_back(SegmentRange(segmentStart, tmp));
            }
            segmentStart = tmp + 1;
            minBackwardResCapacity = backwardResCapacityDP.min_cost_first(tmp + 1, from, tmp);
        }

        if (splitHappened) {
            const auto& oldSegment = backwardResidualSegments.front();
            if (segmentStart < oldSegment.second) {
                splitSegments.push_back(SegmentRange(segmentStart, oldSegment.second));
            }

            // Update segments
            backwardResidualSegments.pop_front();
            backwardResidualSegments.insert(backwardResidualSegments.begin(),
                                            splitSegments.begin(),
                                            splitSegments.end());
        }
    }

    return true;
}

void ForwardBackwardGraph::mergeFlowSolutions(const dp_array<double>& forwardResFlowDP, const dp_array<double>& backwardResFlowDP) {
    std::vector<double> forwardResFlows;
    std::vector<double> backwardResFlows;
    assert(forwardResFlowDP.vectorize(forwardResFlows));
    assert(backwardResFlowDP.vectorize(backwardResFlows));
    for (std::size_t i = 0; i < m_n - 1; ++i) {
        if (forwardResFlows[i] > backwardResFlows[i]) {
            m_forwardEdges[i]->flow = forwardResFlows[i] - backwardResFlows[i];
        } else if (forwardResFlows[i] < backwardResFlows[i]) {
            m_backwardEdges[i]->flow = backwardResFlows[i] - forwardResFlows[i];
        }
    }
}

#ifdef DEBUG_LOTSIZING
void ForwardBackwardGraph::print() const {
    ForwardGraph::print();

    // Inventory backward edges.
    for (const auto& edge: m_backwardEdges) {
        std::cout << "(" << edge->capacity << "," << edge->cost() << "," << edge->flow << ")";
    }
    std::cout << "\n";
}

void ForwardBackwardGraph::print(const std::list<ResidualPath>& forwardResidualPaths,
                                 const std::list<BackwardResidualPathSegment>& backwardResidualPathSegments) const {
    std::cout << "=== Status === \n";
    this->print();
    std::cout << "Forward residual paths:\n";
    ::print(forwardResidualPaths);
    std::cout << "\nBackward residual paths:\n";
    ::print(backwardResidualPathSegments, m_backwardResidualEdgeCostSums);
    std::cout << "\n";
}

void ForwardBackwardGraph::printResiduals() const {
    ForwardGraph::printResiduals();

    // Backward residuals
    for (const auto& edge: m_backwardResidualEdges) {
        std::cout << "(" << edge->capacity << "," << edge->cost() << ")";
    }
    std::cout << "\n";
}
#endif
