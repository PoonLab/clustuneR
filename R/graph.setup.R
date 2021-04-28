#' Makes a graph object based on sequence data and pairwise comparisons
#'
#' Create an implementation of a graph. This is a list, consisting of a distance matrix and some sequence meta data (seq.info)
#' A large part of this process involves resolving growth, ensuring that new sequences are only added prospectively without merging clusters.
#' At this point, we also annotate the minimum retrospective edges of each edge. This is stored in sequence information
#' NOTE: Other growth resolutions may be possible such as random joining, or partial joining, however, these are not currently
#' investigated or implemented
#'
#' @param seq.info: A set of sequence meta-data sorted by alignment header
#' @param edge.info: A pairwise edge matrix of all associated headers in seq.info
#' @param which.new: A set of indices of which sequences were to be labelled "new". this labels certain
#' @param growth.resolution: The method by which growth is resolved. This ensures new cases don't merge clusters
#' By default, each new sequence joins a cluster by only it's minimum retrospective
#' @return A graph, with sequences and edge info. New sequences are only linked by their minimum retrospective edge
#' @export
#' @example examples/create.graph_ex.R
create.graph <- function(seq.info=data.table(), edge.info, which.new=numeric(0), growth.resolution = minimum.retrospective.edge) {

  # Check inputs
  if (nrow(seq.info)==0){
    warning("No sequence meta-data included, creating default seq.info input from headers in edge.info")
    seq.info <- data.table("Header"=colnames(edge.info))
  }
  if (!all(colnames(edge.info) %in% seq.info$Header)) {
    stop("The pairwise distance matrix does not contain all recognized headers from alignment")
  }

  # Assemble graph object
  g <- list()
  g$seq.info <- seq.info
  g$edge.info <- edge.info

  g$seq.info[,"New" := F]
  g$seq.info[which.new,"New" := T]

  g$growth.resolved <- growth.resolution(g)

  return(g)
}

#' A growth resolution helper.
#'
#' Ensures that new sequences only join old clusters through their minimum, retrospective edge
#'
#' @param g: The input graph
#' @return A data table matching each new sequence with its closest retrospective neighbour
minimum.retrospective.edge <- function(g) {

  # Find the minimum retrospective edge of each sequence
  new.seqs <- g$seq.info[(New), Header]
  old.seqs <- g$seq.info[!(New), Header]
  retro.edges <- g$edge.info[new.seqs, old.seqs]
  min.retro.edges <- sapply(seq_len(nrow(retro.edges)), function(i) {
    names(which.min(retro.edges[i, ]))
  })

  if(length(min.retro.edges)==0){
    min.retro.edges <- numeric(0)
  }

  DT <- data.table("NewHeader" = new.seqs, "OldHeader" = min.retro.edges)

  return(DT)
}
