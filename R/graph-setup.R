#' Make a graph object from an edge list and node list.
#'
#' A graph consists of a node list including sequence metadata (`seq.info`) and 
#' an edge list (`edge.info`) of genetic distances between sequences.  
#' clustuneR applies a distance threshold to induce a subgraph comprising one 
#' or more connected components as 'clusters'.
#' 
#' New sequences (identified by `which.new`) are limited to a single in-edge 
#' from one older sequence to ensure that cluster growth does not result in 
#' merging clusters, using the `growth.resolution` method.  Clusters comprising 
#' only new sequences are ignored for subsequent analysis.
#' 
#' @param edge.info: data.frame, containing three columns (ID1, ID2, Distance). 
#'                   ID1 and ID2 must match Header values in seq.info.  TN93 
#'                   (https://github.com/veg/tn93) generates this format by 
#'                   default.
#' @param seq.info: data.table or data.frame, A node list of sequence metadata. 
#'                  Must contain a Header column of sequence labels.
#' @param which.new: numeric or logical, which sequences in seq.info are "new". 
#'                   If numeric, this is a list of integer indices.
#' @param growth.resolution: function, The method by which growth is resolved. 
#'                           This ensures new cases don't merge clusters. By 
#'                           default, each new sequence joins a cluster by only 
#'                           its minimum retrospective edge. 
#' @return S3 object of class clusData
#' @export
read.edges <- function(edge.info, seq.info, which.new=numeric(0), 
                         growth.resolution = minimum.retrospective.edge) {
  # Check inputs
  if (!is.element("Header", names(seq.info)))
    stop("seq.info must contain column `Header`.")
  if (!all(c("ID1", "ID2", "Distance") %in% names(edge.info))) 
    stop("edge.info must contain columns `ID1`, `ID2` and `Distance`.")
  
  edge.labels <- unique(c(edge.info$ID1, edge.info$ID2))
  if (!all(edge.labels %in% seq.info$Header)) {
    stop("edge.info contains ID1 or ID2 strings that are not found in", 
         "seq.info$Header")
  }

  obj <- list()
  obj$seq.info <- as.data.table(seq.info)
  
  # make edge list more compact by replacing names with indices to seq.info
  obj$edge.info <- edge.info
  obj$edge.info$ID1 <- match(obj$edge.info$ID1, seq.info$Header)
  obj$edge.info$ID2 <- match(obj$edge.info$ID2, seq.info$Header)
  
  # exclude self-edges
  obj$edge.info <- obj$edge.info[obj$edge.info$ID1 != obj$edge.info$ID2, ]
  #TODO: remove duplicate edges (a-b and b-a)

  # add a column indicating which sequences are new
  if ("New" %in% names(obj$seq.info))
    warning("Warning, column `New` in seq.info is being overwritten.")
  if (is.numeric(which.new)) {
    obj$seq.info$New <- FALSE
    obj$seq.info$New[which.new] <- TRUE
  } 
  else if (is.logical(which.new)) {
    if (length(which.new) != nrow(obj$seq.info))
      stop("Length mismatch between which.new and seq.info.")
    obj$seq.info$New <- which.new
  } 
  else {
    stop("Unrecognized object type ", typeof(which.new), " for which.new")
  }

  if (!any(obj$seq.info$New))
    stop("At least one sequence must be marked as New, detected none.")
  
  # save resolved single edges between old and new cases
  resolved.edges <- obj$edge.info[growth.resolution(obj), ]
  
  # remove all edges involving new cases
  new.seqs <- which(obj$seq.info$New)
  obj$edge.info <- obj$edge.info[!(obj$edge.info$ID1 %in% new.seqs | 
                                     obj$edge.info$ID2 %in% new.seqs), ]
  
  # append resolved edges
  obj$edge.info <- rbind(obj$edge.info, resolved.edges)
  class(obj) <- "clusData"
  return(obj)
}


#' Generate edge list from a sequence alignment.
#' @param seqs: object of class ape::DNAbin


#' The default growth resolution helper.
#'
#' Ensures that new sequences only join old clusters through the shortest 
#' retrospective edge, i.e., an edge connecting the new node to a non-new node.
#'
#' @param obj: S3 object of class clusData from create.graph
#' @return integer, vector of row indices into g$edge.info
minimum.retrospective.edge <- function(obj) {
  # Find the minimum retrospective edge of each sequence
  new.seqs <- which(obj$seq.info$New)
  old.seqs <- which(!obj$seq.info$New)
  
  # extract row index for shortest edge for each new node
  min.retro.edges <- sapply(new.seqs, function(new.seq) {
    idx <- c(which(obj$edge.info$ID1 == new.seq & obj$edge.info$ID2 %in% old.seqs), 
             which(obj$edge.info$ID2 == new.seq & obj$edge.info$ID1 %in% old.seqs))
    idx[which.min(obj$edge.info$Distance[idx])]
  })
  
  return(min.retro.edges)
}

#TODO: what are some other ways of resolving growth?  random?  closest in date?

print.clusData <- function(obj) {
  cat("clusData S3 object (")
  cat(paste(nrow(obj$seq.info), " nodes, ", nrow(obj$edge.info), " edges)\n\n", sep=""))
  cat(paste(sum(obj$seq.info$New)), " new nodes\n\n", sep="")
  cat("seq.info:\n")
  print(head(obj$seq.info))
  cat("\nedge.info:\n")
  print(head(obj$edge.info))
}
