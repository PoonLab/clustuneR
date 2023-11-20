#' Makes a graph object based on sequence data and pairwise comparisons
#'
#' Create an implementation of a graph. This is a list, consisting of a node
#' list including metadata (seq.info) and an edge list (edge.info). 
#' Any new sequences are resolved as growth ensuring that new sequences are 
#' added prospectively without merging clusters. 
#' Clusters containing purely new sequences are also ignored for this 
#' purpose. At this point, new sequences only connect to the closest non-new 
#' neighbour (although other options for growth resolutions may exist and be 
#' implemented).
#' 
#' @param seq.info: data.table, A set of sequence metadata. Must contain a 
#' Header column of sequence labels.  If none given, then attempts to 
#' extract metadata from headers in edge.info.
#' @param edge.info: data.frame, must contain three columns (ID1, ID2, Distance). 
#' ID1 and ID2 must match Header values in seq.info.  TN93 generates this format 
#' by default.
#' @param which.new: numeric or logical, which sequences in seq.info are "new". 
#' If numeric, this is a list of indices.
#' @param growth.resolution: function, The method by which growth is resolved. This ensures 
#' new cases don't merge clusters. By default, each new sequence joins a cluster 
#' by only it's minimum retrospective edge. 
#' @return S3 object of class clusData
#' See graph.ex for a more detailed example with annotated data.
#' @export
read.edges <- function(seq.info, edge.info, which.new=numeric(0), 
                         growth.resolution = minimum.retrospective.edge) {
  # Check inputs
  if (!all(c("ID1", "ID2", "Distance") %in% names(edge.info))) {
    stop("edge.info must contain columns `ID1`, `ID2` and `Distance`.")
  }
  if (!all(unique(c(edge.info$ID1, edge.info$ID2)) %in% seq.info$Header)) {
    stop("edge.info contains ID1 or ID2 strings that are not found in seq.info$Header")
  }

  obj <- list()
  obj$seq.info <- as.data.table(seq.info)
  
  # make edge list more compact by replacing names with indices to seq.info
  obj$edge.info <- edge.info
  obj$edge.info$ID1 <- match(obj$edge.info$ID1, seq.info$Header)
  obj$edge.info$ID2 <- match(obj$edge.info$ID2, seq.info$Header)

  # add a column indicating which sequences are new
  if (is.numeric(which.new)) {
    obj$seq.info$New <- FALSE
    obj$seq.info$New[which.new] <- TRUE
  } 
  else if (is.logical(which.new)) {
    if (length(which.new) != nrow(obj$seq.info)) {
      stop("Length mismatch between which.new and seq.info.")
    }
    obj$seq.info$New <- which.new
  } 
  else {
    stop("Unrecognized object type ", typeof(which.new), " for which.new")
  }
  
  if (!any(obj$seq.info$New)) {
    stop("At least one sequence must be marked as New, detected none.")
  }
  class(obj) <- "clusData"
  
  # save resolved single edges between old and new cases
  resolved.edges <- obj$edge.info[growth.resolution(obj), ]
  
  # remove all edges involving new cases
  new.seqs <- which(obj$seq.info$New)
  obj$edge.info <- obj$edge.info[!(obj$edge.info$ID1 %in% new.seqs | 
                                     obj$edge.info$ID2 %in% new.seqs), ]
  
  # append resolved edges
  obj$edge.info <- rbind(obj$edge.info, resolved.edges)
  
  return(obj)
}


#' The default growth resolution helper.
#'
#' Ensures that new sequences only join old clusters through the shortest 
#' retrospective edge, i.e., an edge connecting the new node to a non-new node.
#'
#' @param obj: S3 object of class clusData from create.graph
#' @return integer, vector of row indices into g$edge.info
minimum.retrospective.edge <- function(obj) {
  stopifnot(class(obj) == "clusData")

  # Find the minimum retrospective edge of each sequence
  new.seqs <- which(obj$seq.info$New)
  old.seqs <- which(!obj$seq.info$New)
  
  # row indices for edges from an old node to a new node
  retro.edges <- c(which(obj$edge.info$ID1 %in% new.seqs & obj$edge.info$ID2 %in% old.seqs),
                   which(obj$edge.info$ID1 %in% old.seqs & obj$edge.info$ID2 %in% new.seqs))
  
  # extract row index for shortest edge for each new node
  min.retro.edges <- sapply(new.seqs, function(new.seq) {
    my.edges <- c(which(obj$edge.info$ID1 == new.seq | obj$edge.info$ID2 == new.seq))
    my.subset <- obj$edge.info[my.edges, ]
    my.edges[which.min(my.subset$Distance)]
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
