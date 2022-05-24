# clustuneR
### Implementing clustering algorithms on genetic data and finding optimal parameters through the performance of predictive growth models.

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



#### Installation

> Because clustuneR uses [pplacer](https://github.com/matsen/pplacer/) to graft new sequences onto a phylogenetic tree, it can currently only be run on Linux systems.

If you have the [`git`](https://git-scm.com/) version control system installed on your computer, you can clone the repository by navigating to a location of your filesystem where the package will be copied, and then running
```
git clone https://github.com/PoonLab/clustuneR.git 
```

If you do not have `git` installed, then you can download the most recent (developmental version) package as a ZIP archive at this link:
https://github.com/PoonLab/clustuneR/archive/refs/heads/master.zip

or from the Releases page:
https://github.com/PoonLab/clustuneR/releases

If you have downloaded a `.zip` or `.tar.gz` archive, you can use `unzip` or `tar -zvxf` on the command line, or double-click on the archive file in your desktop environment.

Use `cd clustuneR` to enter the package directory and run the following command to install the package into R:
```
R CMD INSTALL .
```
You should see something like this on your console:
```
* installing to library ‘/Library/Frameworks/R.framework/Versions/4.0/Resources/library’
* installing *source* package ‘clustuneR’ ...
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (clustuneR)
```


## Example Usage

#### Setup

To start, we can pull in a tree file and sequence alignment using ape and sequence meta-data as a data-table.
In order to use pplacer and extend the tree, the tree should be built on a subset of the sequences in the alignment, *i.e.*, excluding sequences from the newest year.

```R
setwd("~/git/clustuneR")
require(clustuneR)  # pkgload::load_all()

phy <- ape::read.tree("data/old.treefile")
seqs <- ape::read.FASTA("data/na.fasta", type="DNA")

# parse sequence headers (alternatively import from another file)
seq.info <- pull.headers(seqs, sep="_", var.names=c('accession', 'coldate', 'subtype'),
var.transformations=c(as.character, as.Date, as.factor))

# use pplacer to graft new sequences onto old tree
phy.extend <- extend.tree(phy, seq.info, seqs, mc.cores=4, log.file="data/old.log")

# generate cluster sets under varying parameter settings
param.list <- lapply(seq(0.001, 0.04, 0.001), function(x) list(t=phy.extend, branch.thresh=x, boot.thresh=0.95))

# has variables Size, Growth, 
cluster.sets <- multi.cluster(step.cluster, param.list) 


p.models = list(
    "NullModel" = function(x){
        glm(Growth~Size, data=x, family="poisson")
    },
    "TimeModel" = function(x){
        glm(Growth~Size+coldate, data=x, family="poisson")
    }
)
p.trans = list(
    "coldate" = function(x){mean(x)}
)

res <- fit.analysis(cluster.sets, predictive.models=p.models, 
                    predictor.transformations=p.trans)
AICs <- get.AIC(res)
AIC_loss <- AICs$TimeModelAIC - AICs$NullModelAIC
```

old example:
```
t <- ape::read.tree(<PATH_TO_MY_TREE>)
seqs <- ape::read.FASTA(<PATH_TO_MY_ALIGNMENT>, type = "DNA")
seq.info <- data.table::fread(<PATH_TO_MY_SEQUENCE_METADATA>)

t.extended <- extend.tree(
    t, seq.info, seqs, mc.cores = 4, 
    log.file = "/home/cchato/Data/paperData/tnTrees/tn_Old.fasta.log",
)

# Alternatively - sequence info can be extracted from the headers.
# For headers looking like `UNIQUE-ID_2018_2016_35`
pull.headers(
    seqs, sep="_"
    var.names = c("Sample_ID", "Collection_Year", "Diagnostic_Year", "Age"), 
    var.transformations = list(as.character, as.numeric, as.numeric, as.numeric)
)

param.list <- list(
    list(t = t.extended, branch.thresh = 0.040, boot.thresh = 0.95),
    list(t = t.extended, branch.thresh = 0.030, boot.thresh = 0.95),
    list(t = t.extended, branch.thresh = 0.015, boot.thresh = 0.95),
    list(t = t.extended, branch.thresh = 0.040, boot.thresh = 0),   
    list(t = t.extended, branch.thresh = 0.030, boot.thresh = 0),
    list(t = t.extended, branch.thresh = 0.015, boot.thresh = 0),
)
cluster.sets <- multi.cluster(step.cluster, param.list) 
```

At this point the tree contains

- **t.extended$growth.info** Summaries of pplacer's sequence placements on the tree (branchlengths, confidence, neighbouring nodes).

- **t.extended$node.info** Summaries of node info max and mean patristic distances under a given node of the tree as well as descendant lists.

- **t.extended$seq.info** Sequence metadata

- **t.extended$path.info** Node-paths between nodes. This is used by some clustering functions


#### Clustering

Different clustering methods are under development within clustuneR. For trees using pplacer, the `step.cluster` function should be used.
The parameters are modifiable for each clustering method and a fast `multi.cluster` function can take in many many parameter sets at once to output many sets of clusters.

```R
cluster.set1 <- step.cluster(t = t.extended, branch.thresh = 0.040, boot.thresh = 0.95)
cluster.set2 <- step.cluster(t = t.extended, branch.thresh = 0.030, boot.thresh = 0.95)
cluster.set3 <- step.cluster(t = t.extended, branch.thresh = 0.015, boot.thresh = 0.95)
cluster.set4 <- step.cluster(t = t.extended, branch.thresh = 0.040, boot.thresh = 0)
cluster.set5 <- step.cluster(t = t.extended, branch.thresh = 0.030, boot.thresh = 0)
cluster.set6 <- step.cluster(t = t.extended, branch.thresh = 0.015, boot.thresh = 0)

# Using multi.cluster
param.list <- list(
    list(t = t.extended, branch.thresh = 0.040, boot.thresh = 0.95),
    list(t = t.extended, branch.thresh = 0.030, boot.thresh = 0.95),
    list(t = t.extended, branch.thresh = 0.015, boot.thresh = 0.95),
    list(t = t.extended, branch.thresh = 0.040, boot.thresh = 0),   
    list(t = t.extended, branch.thresh = 0.030, boot.thresh = 0),
    list(t = t.extended, branch.thresh = 0.015, boot.thresh = 0),
)
cluster.sets <- multi.cluster(step.cluster, param.list) 
```

A cluster set will include their size, descendant lists and their Growth (if pplacer was used to propose placements on fixed clusters)

#### GAIC Analysis

Finally, GAIC Analysis can be performed to measure the information content of each cluster set. 
A minimal example involves looking at the AIC loss obtained when switching from a model using only size to a model using size as well as a mean measurement of time.

```R
predictive.models = list(
    "NullModel" = function(x){
        glm(Growth~Size, data=x, family="poisson")
    },
    "TimeModel" = function(x){
        glm(Growth~Size+Diagnostic_Year, data=x, family="poisson")
    }
)
predictor.transformations = list(
    "coldate" = function(x){mean(x)}
)

res <- fit.analysis(cluster.sets, predictive.models, predictor.transformations)
AICs <- get.AIC(res)
AIC_loss <- AICs$TimeModelAIC - AICs$NullModelAIC
```


## References
This package includes the binaries for pplacer and guppy (https://matsen.fhcrc.org/pplacer, released under the GPLv3 license), which are used to add new tips onto a fixed tree to simulate cluster growth prospectively. 

* Matsen FA, Kodner RB, Armbrust EV. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC bioinformatics. 2010 Dec;11(1):1-6.

As an example, this package includes a subset of a larger published HIV-1 *pol* sequence data set. These sequences were originally published in a study by Vrancken *et al.* (2017) and publicly accessible in the GenBank database under the PopSet accession `1033910942`.

* Benson DA, Karsch-Mizrachi I, Lipman DJ, Ostell J, Rapp BA, Wheeler DL. GenBank. Nucleic acids research. 2000 Jan 1;28(1):15-8.

* Vrancken B, Adachi D, Benedet M, Singh A, Read R, Shafran S, Taylor GD, Simmonds K, Sikora C, Lemey P, Charlton CL. The multi-faceted dynamics of HIV-1 transmission in Northern Alberta: A combined analysis of virus genetic and public health data. Infection, Genetics and Evolution. 2017 Aug 1;52:100-5.

