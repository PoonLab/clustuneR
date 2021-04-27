# PhydingClusters (name pending)
Implementing phylogenetic clustering algorithms and finding their optimal parameters through predictive growth modelling.

This package acts on ape's implementation of sequence alignments or phylogenetic trees, sorting taxa into groups based on several different clustering algorithms.
Both graph and tree-based clustering algorithms are implemented as functions, with additional tools to create graphs and extend trees with necessary information. 
The graph and tree setup functions also take header-paired meta-data as input, which can associated with clusters for future analysis.

### Key Features
* Any clustering algorithm reports clusters in the same format (a data.table object). This allows for easy comparison of cluster sets built using different algorithms, or built using the same algorithm with differing parameters.
* Certain sequences may be designated as "New" in order to simulate the growth of clusters. All algorithms implemented in this package have a built-in method of growth resolution, which defines how they assign new cases onto known clusters. Because each cluster can then have an associated growth measurement, predictive models can be fit to cluster-associated meta-data. Comparing the performance of these models on one cluster set vs. another can inform the choice of parameters used to assign cases into clusters.
* Many of these functions take other functions as input (for instance - a specific clustering algorithm is a function and often an input to higher level comparison functions. A particular method of resolving cluster growth may also act as an input for a given clustering algorithm). This design feature is intended to keep the package extensible in the face of an ever increasing number of changes and alterations to genetic clustering algorithms.

### References
This package includes the binaries for pplacer and guppy https://matsen.fhcrc.org/pplacer, which are used to achieve the addition of new tips onto a fixed tree. 
* Matsen FA, Kodner RB, Armbrust EV. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC bioinformatics. 2010 Dec;11(1):1-6.

As an example, this package uses a subset data of a published HIV1 sequence data set. These sequences were originally studied in publication by Vrancken, et al (2017) and accessible publically on genbank under the popset# 1033910942
* Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Rapp BA, Wheeler DL. GenBank. Nucleic acids research. 2000 Jan 1;28(1):15-8.
* Vrancken B, Adachi D, Benedet M, Singh A, Read R, Shafran S, Taylor GD, Simmonds K, Sikora C, Lemey P, Charlton CL. The multi-faceted dynamics of HIV-1 transmission in Northern Alberta: A combined analysis of virus genetic and public health data. Infection, Genetics and Evolution. 2017 Aug 1;52:100-5.
