# Efficient Algorithms for Dynamic Lot-Sizing Problems

This repo implements efficient algorithms to solve the dynamic lot-sizing problem. Dynamic lot-sizing problem is a production planning problem where one needs to decide the production and inventory controls to meet the demands over a period time. For a period of $n$ days, this problem can be graphically illustrated as follows:
<img src="LotSizingProb.PNG" alt="Lot Sizing Problem" width="500"/>

Each arc has a cost $cost$ and a capacity $cap$ properties, where $cost$ is used to compute the cost along a production/inventory-forward/inventory-backward arc, and $cap$ limits the amounts of production/inventory on the respective arcs.

Two cases are considered:
1. All costs are linear functions, where $cost_p$, $cost_{inventory_f}$, and $cost_{inventory_b}$ represent the respective linear coefficients of the arcs.
2. Convex costs, where $cost_p(\cdot)$, $cost_{inventory_f}(\cdot)$, and $cost_{inventory_b}(\cdot)$ represent the convex cost functions on the respective flows along the arcs.
