README
================

# clustuneR

### Optimizating genetic clustering methods on the performance of predictive growth models.

A genetic cluster is a grouping of sequences that are markedly more
similar to each other than to other sequences in the data set. Genetic
clustering has many applications in biology, such as defining taxonomic
groups. In the molecular epidemiology of infectious diseases, it can be
used to characterize the transmission of a pathogen through the
population. For instance, a cluster of sequences can represent a recent
transmission outbreak, especially for rapidly-evolving pathogens.

Most clustering methods require the user to select one or more criteria
defining groups. **clustuneR** provides a statistical framework to
select optimal clustering criteria, based on the premise that the most
effective clustering should maximize our ability to predict the
distribution of new cases among clusters.

### Installation

> ⚠️ Because clustuneR uses
> [pplacer](https://github.com/matsen/pplacer/) to graft new sequences
> onto a phylogenetic tree, it can currently only be run on Linux and
> macOS systems.

#### Download the package

If you have the [`git`](https://git-scm.com/) version control system
installed on your computer, you can clone the repository by navigating
to a location of your filesystem where the package will be copied, and
then running

    git clone https://github.com/PoonLab/clustuneR.git 

If you do not have `git` installed, then you can download the most
recent (developmental version) package as a ZIP archive at this link:
<https://github.com/PoonLab/clustuneR/archive/refs/heads/master.zip>

or from the Releases page:
<https://github.com/PoonLab/clustuneR/releases>

If you have downloaded a `.zip` or `.tar.gz` archive, you can use
`unzip` or `tar -zvxf` on the command line, or double-click on the
archive file in your desktop environment.

#### Running macOS binaries

macOS will prevent you from running the *pplacer* and *guppy* binaries
that are distributed with this package. To bypass this safety mechanism,
you will need to follow these steps:

1.  Use Finder or Terminal to navigate to the `inst` folder in the
    package directory.
2.  Attempt to execute the `pplacer.Darwin` binary. If using Finder,
    double-click on the `pplacer` file, which spawns a Terminal window
    running the binary. If using Terminal, enter the command
    `./pplacer.Darwin`. Your system should display a pop-up with the
    message
    `"pplacer" cannot be opened because the developed cannot be verified`.
    Click on the *Cancel* button to dismiss the pop-up.
3.  Open the System Settings app and click on the Privacy & Security
    tab. In one of the panels, you should see the following label:
    `"pplacer.Darwin" was blocked from use because it is not from an identified developer`.
    Click on the *Allow Anyway* button.
4.  Repeat step 2. The pop-up message should now be changed to
    `macOS cannot verify the developer of "pplacer.Darwin". Are you sure you want to open it?`
    Click on the *Open* button. Your Terminal window should now be
    updated with the following text:  
    `Warning: pplacer couldn't find any sequences to place. Please supply an alignment with sequences to place as an argument at the end of the command line.`  
    This means the program is running properly.
5.  Repeat steps 2-4 by substituting `guppy.Darwin` for
    `pplacer.Darwin`.

You can find similar instructions on the Apple website at
<https://support.apple.com/en-ca/HT202491>

If you do not want to trust the binaries in this package distribution,
you can download the macOS binaries directly from the Matsen lab [GitHub
release
page](https://github.com/matsen/pplacer/releases/tag/v1.1.alpha17) or
compile them from source yourself.

#### Package installation

Use `cd clustuneR` to enter the package directory and run the following
command to install the package into R:

    R CMD INSTALL .

You should see something like this on your console:

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

The process will pause at `moving datasets to lazyloadDB` because there
are several large binary files (`pplacer` and `guppy`) that are included
with this package distribution.

## Usage

**clustuneR** can optimize clustering methods based on either pairwise
genetic distances (graph-based clustering) or a phylogenetic tree
(subtree-based clustering). At minimum, you will need to start with an
alignment of genetic sequences and the respective sample collection
dates as metadata.

The general workflow is:

1.  Partition the sequences into subsets of known and new cases, based
    on collection dates.
2.  Generate sets of clusters from the known cases under varying
    clustering criteria.
3.  Obtain the distribution of new cases among clusters under the
    different criteria as a measure of cluster growth.
4.  Fit a null model of cluster growth as a count outcome predicted by
    cluster size.
5.  Fit an alternate model of cluster growth incorporating additional
    metadata, e.g., sampling dates.
6.  Generate a ∆AIC profile comparing the fits of these models under
    varying clustering criteria.

### Graph-based clustering

A graph consists of a set of nodes and edges. Each node represents an
infection. An edge between nodes indicates that the genetic similarity
of the respective infections falls below some threshold. Conventionally,
we interpret each connected component of the graph as a cluster. A
connected component is a group of nodes such that (1) every node can be
reached from another node through a path of edges, and (1) there are no
edges to nodes outside of the group. Varying the threshold for edges
yields different sets of clusters.

In the following example, we start by reading in a sequence alignment (a
published set of anonymized HIV-1 sequences from Canada), extracting
metadata from the sequence labels, and identifying a subset of new
sequences:

``` r
require(clustuneR)
seqs <- ape::read.FASTA("data/na.fasta", type="DNA")

# parse sequence headers (alternatively import from another file)
seq.info <- parse.headers(names(seqs), sep="_", 
  var.names=c('accession', 'coldate', 'subtype'),
  var.transformations=c(as.character, as.Date, as.factor)
)
seq.info$colyear <- year(seq.info$coldate)

# determine newest year for growth
which.new <- which(seq.info$colyear == max(seq.info$colyear))
```

> You may already have these metadata in the form of a tabular data set
> (*i.e.*, a CSV file), in which case you can simply load these metadata
> as a data frame.

Next, we need to load a list of edges, where each row specifies two node
labels and a distance. These data can be generated from a sequence
alignment using the program [TN93](https://github.com/veg/tn93). The
resulting output file is enormous (\>34MB), so we do not include it in
this package!

``` r
# load genetic distances (run `tn93 -t 1 -o na.tn93.csv na.fasta`)
edge.info <- read.csv("data/na.tn93.csv")
obj <- read.edges(edge.info, seq.info, which.new)

# generate cluster sets under varying parameter settings
cutoffs <- seq(0, 0.04, length.out=50)
param.list <- lapply(cutoffs, function(x) { 
  list(dist.thresh=x, time.var="colyear") 
  })
cluster.sets <- multi.cluster(obj, param.list, component.cluster) 
```

By specifying a `time.var` argument in `param.list`, we are fitting a
model to the distribution of sample collection years to predict edges
between cases. For a more detailed explanation of this method, please
refer to the vignettes.

The last step of the analysis is to fit regression models to the
distribution of new cases among clusters.

``` r
ptrans <- list("Weight"=sum)
pmods <- list(
  "NullModel"=function(x) glm(Growth~Size, data=x, family="poisson"),
  "AltModel"=function(x) glm(Growth~Weight, data=x, family="poisson")
)
res <- fit.analysis(cluster.sets, models=pmods, transforms=ptrans)
gaic <- get.AIC(res, param.list)
```

Here, `gaic` is a data frame that stores the key result of our
analysis - the AIC values associated with the two models under varying
clustering thresholds. The optimal TN93 distance cutoff is identified by
the greatest difference between the AICs of the alternative and null
models, which we can visualize as a plot:

``` r
par(mar=c(5,5,1,1))
plot(cutoffs, gaic$AltModel - gaic$NullModel, type='l', 
     lwd=2, col='cadetblue',
     xlab="TN93 distance cutoffs", ylab="delta-AIC")
abline(h=0, lty=2)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Tree-based clustering

A phylogenetic tree is a hypothesis about how different populations are
related by their common ancestors. In the context of molecular
epidemiology, the ancestral nodes in a tree relating different
infections can approximate transmission events in the past. Thus, a
cluster of sequences connected by short branches in the tree may
represent an outbreak.

As above, we start the same set of anonymized HIV-1 sequences in
`data/na.fasta`. First, we generate a new alignment excluding any
sequences collected in the most recent year:

``` r
ape::write.FASTA(seqs[-which.new], file="data/na-old.fasta")
```

Next, we use a maximum likelihood program such as
[IQ-TREE](http://www.iqtree.org/) to reconstruct a tree relating these
“old” sequences:

``` console
iqtree -bb 1000 -m GTR -nstop 200 -s na-old.fasta
```

Note we’ve requested a specific model of nucleotide substitution (GTR)
to bypass the model selection stage of this program. Even so, this is a
time-consuming step - to speed things up, we’ve provided these IQ-TREE
output files at `data/na.nwk` and `data/na.log`.

> **clustuneR** uses a program (`pplacer`) that can work with the
> outputs of IQ-TREE,
> [FastTree](http://www.microbesonline.org/fasttree/) and
> [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/). You’ll
> have to specify which ML tree reconstruction program you used in the
> next step.

Assuming you’ve kept R running, our next step is to import the ML tree
into R. (If you quit R, you’ll have to repeat the previous steps to
import the alignment and parse headers.) We can then use `pplacer` to
use maximum likelihood to graft the “new” sequences onto this tree.

``` r
phy <- ape::read.tree("data/na.nwk")
phy <- import.tree(phy, seq.info)
phy.extend <- extend.tree(phy, seqs, log.file="data/na.log")
```

We can reuse the `cutoffs` vector from the previous example to configure
a new parameter list for generating different sets of clusters. In this
case, we have two criteria: (1) a threshold for the total branch length
from each tip to the root of a subtree, and (2) the bootstrap support
for the subtree:

``` r
param.list <- lapply(cutoffs, function(x) list(branch.thresh=x, boot.thresh=0.95))
cluster.sets <- multi.cluster(phy.extend, param.list, step.cluster) 
```

We also need to specify two different regression models to fit to these
sets of clusters. Unlike our graph clustering example, we are going to
simply add the mean sample collection date of sequences in each cluster
as a second model term:

``` r
p.models = list(
  "NullModel"=function(x) glm(Growth~Size, data=x, family="poisson"),
  "TimeModel"=function(x) glm(Growth~Size+coldate, data=x, family="poisson")
)
# average sample collection dates across nodes in each cluster
p.trans = list("coldate"=mean)
res <- fit.analysis(cluster.sets, models=p.models, transforms=p.trans)
gaic <- get.AIC(res, param.list)
```

Finally, we can plot the difference in AIC between the models to select
an optimal branch threshold:

``` r
par(mar=c(5,5,1,1))
plot(cutoffs, gaic$TimeModel - gaic$NullModel, type='l', 
     lwd=2, col='cadetblue', xlab="Branch threshold", ylab="delta-AIC")
abline(h=0, lty=2)
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## References

If you use **clustuneR** for your work, please cite one of the following
references:

- Chato C, Kalish ML, Poon AF. Public health in genetic spaces: a
  statistical framework to optimize cluster-based outbreak detection.
  Virus evolution. 2020 Jan;6(1):veaa011.

- Chato C, Feng Y, Ruan Y, Xing H, Herbeck J, Kalish M, Poon AF.
  Optimized phylogenetic clustering of HIV-1 sequence data for public
  health applications. PLOS Computational Biology. 2022 Nov
  30;18(11):e1010745.

This package includes the binaries for pplacer and guppy
(<https://matsen.fhcrc.org/pplacer>, released under the GPLv3 license),
which are used to add new tips onto a fixed tree to simulate cluster
growth prospectively.

- Matsen FA, Kodner RB, Armbrust EV. pplacer: linear time
  maximum-likelihood and Bayesian phylogenetic placement of sequences
  onto a fixed reference tree. BMC bioinformatics. 2010 Dec;11(1):1-6.

This package includes some anonymized HIV-1 sequences that were placed
in the public domain in association with the following publication:

- Vrancken B, Adachi D, Benedet M, Singh A, Read R, Shafran S, Taylor
  GD, Simmonds K, Sikora C, Lemey P, Charlton CL. The multi-faceted
  dynamics of HIV-1 transmission in Northern Alberta: A combined
  analysis of virus genetic and public health data. Infection, Genetics
  and Evolution. 2017 Aug 1;52:100-5.
