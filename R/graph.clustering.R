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
#' @param time.var:  character, column name for discrete time variable to fit 
#'                   a model of edge density decay with time (optional)
#' @return data.frame, known cases annotated with cluster ID and growth
#' @export
#' @example examples/component.cluster_ex.R
component.cluster <- function(obj, dist.thresh=0, setID=0, time.var=NA) {

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
  cluster.set <- obj$seq.info[
    # carry over columns from seq.info
    !(New), lapply(seq.cols, function(nm) list(get(nm))), by=Cluster
    ]
  cluster.set[, "Size" := length(V1[[1]]), by = 1:nrow(cluster.set)]
  colnames(cluster.set) <- c("ClusterID", seq.cols, "Size")
  cluster.set$New <- NULL  # should be all FALSE
  cluster.set <- cluster.set[order(ClusterID),]
  
  # fit edge probability decay model (e.g., fit.decay="colyear")
  if (!is.na(time.var)) {
    times <- obj$seq.info[[time.var]]
    time.counts <- table(times)
    if (sum(time.counts == 1) / length(time.counts) > 0.5) {
      stop("fit.decay only supports discrete time, ", time.var, 
           " appears to be continuous")
    }
    
    # remove new cases and associated edges
    max.time <- max(times)
    keep <- which(times < max.time)
    old.nodes <- obj$seq.info[keep, ]
    old.nodes$time <- times[keep]
    old.edges <- filtered.edges[filtered.edges$ID1 %in% keep & filtered.edges$ID2 %in% keep, ]
    old.edges$t1 <- times[old.edges$ID1]
    old.edges$t2 <- times[old.edges$ID2]
    
    # fit binomial model to bipartite graph
    fit <- fit.decay(old.nodes, old.edges)
    
    # the predicted edge probability for a cluster is a function of mean node age
    #cluster.set[, "EdgeProb" := ]
  }
  
  # Attach growth info and set ID
  growth <- table(obj$seq.info[(New) & (Cluster %in% cluster.set$ClusterID), Cluster])
  cluster.set[, "Growth" := 0]
  cluster.set[ClusterID %in% as.numeric(names(growth)), Growth := as.numeric(growth)]
  
  cluster.set[, "DistThresh" := dist.thresh]
  cluster.set[, "SetID" := setID]
  
  class(cluster.set) <- c("cluster.set", class(cluster.set))
  
  return(cluster.set)
}


#' @param nodes:  data.table, excluding new cases
#' @param edges:  data.frame, edges below a given distance threshold
fit.decay <- function(nodes, edges) {
  # for every node, find the shortest edge from an older node
  # FIXME: this could be done once only outside this function
  positives <- lapply(1:nrow(nodes), function(child) {
    child.t <- nodes$time[child]
    my.edges <- edges[(edges$ID1==child | edges$ID2==child), ]
    parents <- ifelse(my.edges$ID1==child, my.edges$ID2, my.edges$ID1)
    parents.t <- nodes$time[parents]
    
    my.edges$dt <- child.t - parents.t
    my.edges <- my.edges[parents.t < child.t, ]
    my.edges$dt[which.min(my.edges$Distance)]
  })
  positives <- unlist(positives)
  counts <- as.data.frame(table(positives))
  names(counts) <- c("dt", "positives")
  
  # calculate negative counts
  time.vals <- sort(unique(times), decreasing=TRUE)
  counts$total <- 0
  for (i in 1:(length(time.vals)-1)) {
    t2 <- time.vals[i]
    for (j in (i+1):length(time.vals)) {
      t1 <- time.vals[j]
      dt <- t2 - t1
      possible.edges <- as.integer(
        time.counts[as.character(t1)] * time.counts[as.character(t2)]
      )
      counts$total[counts$dt == dt] <- counts$total[counts$dt == dt] + 
        possible.edges
    }
  }
  
  fit <- glm(cbind(positives, total) ~ as.integer(dt), family="binomial", 
             data=counts)
  return(fit)
}
