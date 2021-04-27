# PhydingClusters (name pending)
Implementing phylogenetic clustering algorithms and finding their optimal parameters through predictive growth modelling.

This package acts on ape's implementation of sequence alignments or phylogenetic trees, sorting taxa into groups based on several different clustering algorithms.
Both graph and tree-based clustering algorithms are implemented as functions, with additional tools to create graphs and extend trees with necessary information. 
The graph and tree setup functions also take header-paired meta-data as input, which can associated with clusters for future analysis.

###Key Features
* Any clustering algorithm reports clusters in the same format (a data.table object). This allows for easy comparison of cluster sets built using different algorithms, or built using the same algorithm with differing parameters.
* Certain sequences may be designated as "New" in order to simulate the growth of clusters. All algorithms implemented in this package have a built-in method of growth resolution, which defines how they assign new cases onto known clusters. Because each cluster can then have an associated growth measurement, predictive models can be fit to cluster-associated meta-data. Comparing the performance of these models on one cluster set vs. another can inform the choice of parameters used to assign cases into clusters.
* Many of these functions take other functions as input (for instance - a specific clustering algorithm is a function and often an input to higher level comparison functions. A particular method of resolving cluster growth may also act as an input for a given clustering algorithm). This design feature is intended to keep the package extensible in the face of an ever increasing number of changes and alterations to genetic clustering algorithms.
