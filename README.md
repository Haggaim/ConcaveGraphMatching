# (Probably) Concave Graph Matching

-----------------
## Abstract 
In this paper we address the graph matching problem. Following the recent works of zaslavskiy2009path,Vestner2017 we analyze and generalize the idea of concave relaxations. We introduce the concepts of conditionally concave and probably conditionally concave energies on polytopes and show that they encapsulate many instances of the graph matching problem, including matching Euclidean graphs and graphs on surfaces. We further prove that local minima of probably conditionally concave energies on general matching polytopes (e.g., doubly stochastic) are with high probability extreme points of the matching polytope (e.g., permutations).

-----------------
## Code for the Neurips 2018 paper 
This code implements graph matching with the Frank-Wolfe algorithm, as described in the paper. The examples match shape from the SHREC07 dataset.
Run matchShapes.m or matchShapesOneSided.m for optimizing over permutations or one-sided permutations repectively. 
