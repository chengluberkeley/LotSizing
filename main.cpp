//
//  main.cpp
//  LotSizing
//
//  Created by ChengLu on 11/24/21.
//  Copyright © 2021 Cheng Lu. All rights reserved.
//

#include "LotSizingSolver.hpp"

#include <chrono>
#include <iostream>
#include <random>

int main(int argc, const char * argv[]) {
    // MARK: - Test forward graphs

    std::cout << "Generate forward graph and compute: " << std::endl;

    // Generate forward graphs with random generation.
    std::size_t n = 100;
    std::size_t iterations = 100;

    // Sample demands
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> demandDist(10, 100);
    std::uniform_int_distribution<uint32_t> deltaDist(0, 10);
    std::uniform_real_distribution<double> fCostDist(10.0, 20.0);
    std::uniform_real_distribution<double> bCostDist(0.1, 0.5);

    for (std::size_t i = 0; i < iterations; ++i) {
        Demands demands(n);
        for (std::size_t i = 0; i < n; ++i) {
            demands[i] = demandDist(gen);
        }

        std::vector<uint32_t> productionCapacities(n);
        std::vector<double> productionCosts(n);
        for (std::size_t i = 0; i < n; ++i) {
            productionCapacities[i] = demands[i] + deltaDist(gen);
            productionCosts[i] = fCostDist(gen);
        }

        std::vector<uint32_t> forwardCapacities(n - 1);
        std::vector<double> forwardCosts(n - 1);
        for (std::size_t i = 0; i < n - 1; ++i) {
            forwardCapacities[i] = demandDist(gen);
            forwardCosts[i] = fCostDist(gen);
        }

        // MARK: - Forward graph

        std::cout << "Finish forward graph instance generation." << std::endl;

        ForwardGraph forwardGraph(n, demands, productionCapacities, productionCosts, forwardCapacities, forwardCosts);
        ForwardGraph fasterForwardGraph(n, demands, productionCapacities, productionCosts, forwardCapacities, forwardCosts);

        auto start = std::chrono::steady_clock::now();
        assert(forwardGraph.solve());
        auto end = std::chrono::steady_clock::now();
        auto elapsedSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Runtime = " << elapsedSeconds.count() << " ms\n";

        double cost0 = forwardGraph.cost();
        std::cout << "The total cost is " << cost0 << std::endl;

        assert(forwardGraph.constraintsSatisfied());
        assert(forwardGraph.isOptimal());

        // MARK: - Forward graph faster version

        start = std::chrono::steady_clock::now();
        assert(fasterForwardGraph.fasterSolve());
        end = std::chrono::steady_clock::now();
        elapsedSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Runtime = " << elapsedSeconds.count() << " ms\n";

        double cost0_1 = fasterForwardGraph.cost();
        std::cout << "The total cost of the faster version is " << fasterForwardGraph.cost() << std::endl;

        assert(fasterForwardGraph.constraintsSatisfied());
        assert(fasterForwardGraph.isOptimal());

        assert(fabs(cost0 - cost0_1) < 1e-6);

        // MARK: - Forward-backward graph

        std::vector<uint32_t> backwardCapacities(n - 1);
        std::vector<double> backwardCosts(n - 1);
        for (std::size_t i = 0; i < n - 1; ++i) {
            backwardCapacities[i] = demandDist(gen);
            backwardCosts[i] = bCostDist(gen);
        }

        std::cout << "Finish forward-backward graph instance generation." << std::endl;

        ForwardBackwardGraph forwardBackwardGraph(n, demands, productionCapacities, productionCosts,
                                                  forwardCapacities, forwardCosts,
                                                  backwardCapacities, backwardCosts);
        ForwardBackwardGraph fasterForwardBackwardGraph(n, demands, productionCapacities, productionCosts,
                                                        forwardCapacities, forwardCosts,
                                                        backwardCapacities, backwardCosts);

        start = std::chrono::steady_clock::now();
        assert(forwardBackwardGraph.solve());
        end = std::chrono::steady_clock::now();
        elapsedSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Runtime = " << elapsedSeconds.count() << " ms\n";

        double cost1 = forwardBackwardGraph.cost();
        std::cout << "The total cost is " << cost1 << std::endl;

        start = std::chrono::steady_clock::now();
        assert(fasterForwardBackwardGraph.fasterSolve());
        end = std::chrono::steady_clock::now();
        elapsedSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Faster runtime = " << elapsedSeconds.count() << " ms\n";

        double cost1_1 = fasterForwardBackwardGraph.cost();
        std::cout << "The total cost of the faster version is " << cost1_1 << std::endl;

        assert(forwardBackwardGraph.constraintsSatisfied());
        assert(forwardBackwardGraph.isOptimal());
        assert(fasterForwardBackwardGraph.constraintsSatisfied());
        assert(fasterForwardBackwardGraph.isOptimal());
        assert(fabs(cost1 - cost1_1) < 1e-6);

        assert(cost1 <= cost0);
    }

    return 0;
}
