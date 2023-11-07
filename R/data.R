#' An alignment of HIV1, subtype B sequences
#'
#' A dataset containing 10 HIV1, subtype B polymerase sequences collected in Northern 
#' Alberta Canada. This a 10 sequence sample from popset# 1033910942 on NCBI's genbank 
#' Archive. The sequence headers also include meta-data for sequences.
#'
#' @format An ape DNA object: 10 DNA sequences in binary format stored in a list. 
#' All sequences of same length: 1017 
#' @source \url{ https://www.ncbi.nlm.nih.gov/popset?DbFrom=nuccore&Cmd=Link&LinkName=nuccore_popset&IdsFromResult=1033912042 }
"alignment.ex"

#' An example set of sequence meta.data corresponding to alignment.ex
#'
#' Built from alignment.ex, the example 10 sequence alignment using parse.headers. 
#' The date of each sequence's collection, it's genbank unique accession ID, and 
#' sequence subtype are referenced within the header
#'
#' @format A data.table object with 4 variables:
#' \describe{
#' \item{ID}{Accession IDs (characters) of sequences}
#' \item{CollectionDate}{Collection date of sequences. Full dates given as yyyy-mm-dd}
#' \item{Subtype}{Subtypes (factors) of sequences}
#' \item{Header}{The original headers from the alignement. This matches meta 
#' data to sequences in original alignment}
#' }
"seq.info.ex"

#' An example set of clusters, built using component.cluster
#'
#' A dataset describing 5 different clusters. The headers (from alignment.ex), and 
#' associated meta data (from seq.info.ex) of cluster members is captured, as well
#' as several cluster-level traits, such as growth and size. See component.cluster 
#' or further information onhow these were assigned based on graph.ex as an input
#'
#' @format A data.table object with 9 variables:
#' \describe{
#' \item{ClusterID}{ The unique identifier number for this cluster. A numeric}
#' \item{ID}{A list of vectors, each containing the accession IDs (characters) 
#' of sequences within a cluster}
#' \item{CollectionDate}{A list of vectors, each containing the collection date 
#' of sequences within a cluster}
#' \item{Subtype}{A list of vectors, each containing the subtypes (factors) within 
#' a cluster}
#' \item{Header}{A list of vectors, each containing the original headers from the 
#' alignement used to build this set of clusters}
#' \item{Size}{The original size of this cluster before being updated with new cases. 
#' This simply the number of sequences within the cluster}
#' \item{Growth}{The growth of the cluster after new cases are added}
#' \item{DistThresh}{The pairwise distance threshold used to create this complete 
#' set of clusters. Corresponds to a setID as an input parameter}
#' \item{SetID}{The unique identifier for this set of clusters. A numeric}
#' }
"cluster.ex"

#' An example graph, built based on pairwise TN93 distances
#'
#' This implementation of a graph is a list, describing a set of sequences and the 
#' distances between them. See create.graph for more information on how this graph 
#' was created using alignment.ex as input
#'
#' @format A list of 3 variables
#' \describe{
#' \item{seq.info}{ See seq.info.ex, a data.table containing sequence meta data}
#' \item{edge.info}{ A named, weighted pairwise distance matrix. }
#' \item{growth.resolved}{ a data.table pairing new sequences, to a single sequence. 
#' In order to ensure that clusters do not merge upon growth, a new sequence may only
#' join up to one old cluster. By default, new sequences join the cluster of the 
#' sequence they are most similar to. }
#' }
"graph.ex"

#' A tree built based on a subset of alignment.ex
#'
#' A maximum likelihood tree built using IQ-TREE with model selection and 1000 
#' parametric bootstraps. The log information for this tree is stored in data/IQTREE_log_ex.txt. 
#' A subset of six older sequences (collected before January 1st 2012) from alignment.ex 
#' was used to construct this tree.
#'
#'
#' @format An unrooted, phylogenetic tree with 6 tips and 4 internal nodes. 
#' Node labels represent certainty. See ape's implementation of phylogenetic tree 
#' objects for information about tags within this object
"old.tree.ex"

#' A tree built from alignment.ex
#'
#' This is a maximum likelihood tree built using IQ-TREE with automatic model 
#' selection and 1000 parametric bootstraps. Contrasting old.tree.ex. This is a 
#' complete tree containing all sequences in alignment.ex
#'
#' @format An unrooted, phylogenetic tree with 10 tips and 8 internal nodes. Node 
#' labels represent certainty. See ape's implementation of phylogenetic tree objects 
#' for information about tags within this object.
"full.tree.ex"

#' An extension of an ape tree object which can be used to create clusters
#'
#' An extension of old.tree.ex maximum likelihood tree built using IQ-TREE with automatic 
#' model selection and 1000 parametric bootstraps. growth information and additional 
#' information useful for cluster identification were added by extend.tree. 
#'
#' @format A phylogenetic tree with 6 tips and 4 internal nodes. Node labels represent 
#' certainty. See ape's implementation of phylogenetic tree objects for information 
#' about tags within this object. In addition, there are 4 new objects created by 
#' functions within tree.setup.R
#' \describe{
#' \item{seq.info}{ See seq.info.ex, a data.table containing sequence meta data}
#' \item{node.info}{ Grouping of the meta.data present in seq.info assigned to 
#' various nodes in the tree, coupled with information important to clustering, 
#' such as mean divergence from root, or node certainty }
#' \item{path.info}{ Information regarding the path of edges from tips to the root 
#' of the tree. This is also necessary for some clustering algorithms, specifically 
#' step.cluster}
#' \item{growth.info}{ a data.table pairing new sequences, to a single node in the 
#' tree based on placements assigned by guppy and pplacer. The certainty of this placement, 
#' terminal branch length, neighbour, and branch length from new internal node to 
#' new neighbour are described}
#' }
"extended.tree.ex"