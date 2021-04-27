#'An alignment of HIV1, subtype B sequences
#'
#'A dataset containing 10 HIV1, subtype B polymerase sequences collected in Northern Alberta Canada.
#'This a 10 sequence sample from popset# 1033910942 on NCBI's genbank Archive
#'
#'@format An ape DNA object: 10 DNA sequences in binary format stored in a list. All sequences of same length: 1017 
#'
#'@source \url{https://www.ncbi.nlm.nih.gov/popset?DbFrom=nuccore&Cmd=Link&LinkName=nuccore_popset&IdsFromResult=1033912042}
"alignment.ex"

#'An example set of clusters, built using component.cluster
#'
#'A dataset describing 5 different clusters. Their member headers are listed, as well as the growth they experienced 
#'(ie. the number of new sequences forming clusters with old sequences.). See component.cluster for further information on 
#'how these were assigned based on graph.ex as an input
#'
#' @format A data.table object with 9 variables:
#' \describe{
#' \item{ClusterID}{ The unique identifier number for this cluster. A numberic}
#' \item{ID} {A list of vectors, each containing the accession IDs (characters) of sequences within a cluster}
#' \item{CollectionDate} {A list of vectors, each containing the collection date of sequences within a cluster}
#' \item{Subtype} {A list of vectors, each containing the subtypes (factors) within a cluster}
#' \item{Header}{A list of vectors, each containing the original headers from the alignement used to build this set of clusters}
#' \item{Size} {The original size of this cluster before being updated with new cases. This simply the number of sequences within the cluster}
#' \item{Growth} {The growth of the cluster after new cases are added}
#' \item{DistThresh} {The pairwise distance threshold used to create this complete set of clusters. Corresponds to a setID as an input parameter}
#' \item{SetID} {The unique identifier for this set of clusters. A numeric}
#' \
#' }
"cluster.ex"

#'An example graph, built based on pairwise TN93 distances
#'
#'This implementation of a graph is a list, describing a set of sequences and the distances between them. 
#'See create.graph for more information on how this graph was created using alignment.ex as input
#'
#'
#' @format A list of 3 variables
#' \describe{
#' \item{seq.info}{ See seq.info.ex, a data.table containing sequence meta data}
#' \item{edge.info}{ A named, weighted pairwise distance matrix. }
#' \item{growth.resolved}{ a data.table pairing new sequences, to a single sequence. 
#' In order to ensure that clusters do not merge upon growth, a new sequence may only
#' join up to one old cluster. By default, new sequences join the cluster of the 
#' sequence they are most similar to. }
#' \
#' }
"graph.ex"