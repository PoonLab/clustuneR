require(igraph)

#' Create clusters based on the components of a graph
#'
#' Edges are filtered away using a distance threshold to break up the completely 
#' connected graph such that only similar edges remain.
#'
#' @param obj: S3 object of class clusData.  Contains vertex, edge and growth 
#' resolution information.
#' @param dist.thresh: double, the maximum distance defining which edges are 
#' filtered. A higher distance threshold implies a larger average cluster size
#' @param setID: A numeric identifier for this cluster set.
#' @return data.frame, known cases annotated with cluster ID and growth
#' @export
#' @example examples/component.cluster_ex.R
component.cluster <- function(obj, dist.thresh=0, setID=0) {

  # Filter edges above the distance threshold and prepare for component finding algorithm
  # All edges from a new sequence are filtered except for their "growth-resolved" edge
  filtered.edges <- obj$edge.info[obj$edge.info$Distance <= dist.thresh, ]
  
  g <- graph_from_edgelist(as.matrix(filtered.edges[c("ID1", "ID2")]), directed=FALSE)
  
  # graph_from_edgelist interprets edge list as vertex IDs so it will tend to 
  # omit some number of nodes with no edges past filter
  if (length(V(g)) < nrow(obj$seq.info)) {
    orphans <- seq(max(V(g))+1, nrow(obj$seq.info))
    g <- add_vertices(g, length(orphans))  
  }
  
  # extract network components
  comps <- components(g)
  
  # label sequences with cluster indices
  obj$seq.info$Cluster <- comps$membership
  
  # group known cases by cluster membership
  seq.cols <- colnames(obj$seq.info)
  cluster.set <- obj$seq.info[!(New), lapply(seq.cols, 
                                             function(nm) list(get(nm))), by=Cluster]
  cluster.set[, "Size" := length(V1[[1]]), by = 1:nrow(cluster.set)]
  colnames(cluster.set) <- c("ClusterID", seq.cols, "Size")
  cluster.set$New <- NULL  # should be all FALSE
  cluster.set <- cluster.set[order(ClusterID),]
  
  # Attach growth info and set ID
  growth <- table(obj$seq.info[(New) & (Cluster %in% cluster.set$ClusterID), Cluster])
  cluster.set[, "Growth" := 0]
  cluster.set[ClusterID %in% as.numeric(names(growth)), Growth := as.numeric(growth)]
  
  cluster.set[, "DistThresh" := dist.thresh]
  cluster.set[, "SetID" := setID]
  
  return(cluster.set)
}
