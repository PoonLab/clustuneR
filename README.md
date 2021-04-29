# clustuneR
Implementing clustering algorithms on genetic data and finding optimal parameters through the performance of predictive growth models.

clustuneR builds clusters from inputted sequence alignments and/or phylogenetic trees, allowing users to choose between multiple cluster-building algorithms implemented in the package.
These algorithms can be further augmented through the selection of parameters, such as a required similarity for cluster formation, or a required level of certainty.
The package also takes in meta-data associated with sequences such as a known collection date or subtype/variant classification.
These can also allow users to identify cluster-level characteristics, such as the range of collection dates or the most common subtype/variant within a cluster.

If a subset of sequences are specified as "New", then clustuneR simulates cluster growth by building clusters in two stages: 
first clusters are built from sequences which are not specified as new, then the new sequences are added to clusters. 
Depending on the clustering method used, this second step may include compromises to insure that new sequences do not retroactively change the membership of clusters. 
For example, if a single new sequence forms a cluster with two, previously separate clusters, then those two clusters would have ambiguous growth.
Pairing cluster-level meta-data, with the growth of clusters is a common goal in research and clustuneR contains some functions to help test predictive models based on cluster data.
Furthermore, clustuneR facilitates the assignment of multiple cluster sets from the same data using different methods and parameters.
Pairing these with the effectiveness of growth models can be useful in method/parameter selection.

### References
This package includes the binaries for pplacer and guppy https://matsen.fhcrc.org/pplacer, which are used to add new tips onto a fixed tree to simulate cluster growth prospectively. 
* Matsen FA, Kodner RB, Armbrust EV. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC bioinformatics. 2010 Dec;11(1):1-6.

As an example, this package includes a subset of a larger published HIV1 pol sequence data set. These sequences were originally studied in publication by Vrancken, et al (2017) and accessible publically on genbank under the popset# 1033910942
* Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Rapp BA, Wheeler DL. GenBank. Nucleic acids research. 2000 Jan 1;28(1):15-8.
* Vrancken B, Adachi D, Benedet M, Singh A, Read R, Shafran S, Taylor GD, Simmonds K, Sikora C, Lemey P, Charlton CL. The multi-faceted dynamics of HIV-1 transmission in Northern Alberta: A combined analysis of virus genetic and public health data. Infection, Genetics and Evolution. 2017 Aug 1;52:100-5.
