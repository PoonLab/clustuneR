#' Create clusters based on the components of a graph
#'
#' Edges are filtered away using a distance threshold to break up the completely 
#' connected graph such that only similar edges remain.
#'
#' @param g: The input graph, annotated with vertex, edge, and growth resolution
#' information
#' @param dist.thresh: The maximum distance defining which edges are filtered.
#' A higher distance threshold implies a larger average cluster size
#' @param setID: A numeric identifier for this cluster set.
#' @return A set of clusters as a data.table. See example cluster.ex object 
#' documentation for an example of clustered sequence data + meta data
#' @export
#' @example examples/component.cluster_ex.R
component.cluster <- function(g, dist.thresh = 0, setID = 0) {

  # Filter edges above the distance threshold and prepare for component finding algorithm
  # All edges from a new sequence are filtered except for their "growth-resolved" edge
  filtered.edges <- g$edge.info <= dist.thresh
  diag(filtered.edges) <- F
  filtered.edges[which(g$seq.info$New), ] <- F
  filtered.edges[,which(g$seq.info$New)]  <- F
  filtered.edges[g$growth.resolved$NewHeader, g$growth.resolved$OldHeader] <-
    g$edge.info[g$growth.resolved$NewHeader, g$growth.resolved$OldHeader] <= dist.thresh
  filtered.edges[g$growth.resolved$OldHeader, g$growth.resolved$NewHeader] <-
    g$edge.info[g$growth.resolved$OldHeader, g$growth.resolved$NewHeader] <= dist.thresh
  
  # Run homogenization algorithm to label sequences with their cluster
  seq.cols <- colnames(g$seq.info)
  previous.cluster <- rep(0, nrow(g$seq.info))
  g$seq.info[, "Cluster" := 1:nrow(g$seq.info)]

  while (any(g$seq.info$Cluster != previous.cluster)) {
    previous.cluster <- g$seq.info$Cluster
    g$seq.info[, Cluster := sapply(1:nrow(g$seq.info), function(i) {
      x <- g$seq.info[which(filtered.edges[i, ]), Cluster]
      if (length(x) == 0) {
        return(i)
      } else {
        return(min(x))
      }
    })]
    print(g$seq.info[which(g$seq.info$Cluster != previous.cluster)])
  }

  cluster.set <- g$seq.info[!(New), lapply(seq.cols, function(nm) {
    list(get(nm))
  }), by = Cluster]
  cluster.set[, "Size" := length(V1[[1]]), by = 1:nrow(cluster.set)]
  colnames(cluster.set) <- c("ClusterID", seq.cols, "Size")
  cluster.set$New <- NULL
  
  cluster.set <- cluster.set[order(ClusterID),]
  
  # Attach growth info and set ID
  growth <- table(g$seq.info[(New) & (Cluster %in% cluster.set$ClusterID), Cluster])
  cluster.set[, "Growth" := 0]
  cluster.set[ClusterID %in% as.numeric(names(growth)), Growth := as.numeric(growth)]

  cluster.set[, "DistThresh" := dist.thresh]
  cluster.set[, "SetID" := setID]

  return(cluster.set)
}
