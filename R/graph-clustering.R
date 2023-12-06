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
  # Filter edges above the distance threshold 
  filtered.edges <- obj$edge.info[obj$edge.info$Distance <= dist.thresh, ]
  
  # unconnected vertices will be induced by maximum numeric vertex ID of edgelist
  g <- graph_from_edgelist(as.matrix(filtered.edges[c("ID1", "ID2")]), 
                           directed=FALSE)
  
  # append vertices with numeric IDs above maximum ID in edgelist
  if (length(V(g)) < nrow(obj$seq.info)) {
    orphans <- seq(max(V(g))+1, nrow(obj$seq.info))
    g <- add_vertices(g, length(orphans))  
  }
  
  # extract connected components from graph
  comps <- components(g)
  
  # label sequences with cluster indices
  obj$seq.info$Cluster <- comps$membership
  
  # carry over columns from seq.info
  seq.cols <- colnames(obj$seq.info)
  cluster.set <- obj$seq.info[
    !(New), lapply(seq.cols, function(nm) list(get(nm))), by=Cluster
    ]
  
  # group known cases by cluster membership
  cluster.set[, "Size" := length(V1[[1]]), by = 1:nrow(cluster.set)]
  colnames(cluster.set) <- c("ClusterID", seq.cols, "Size")
  cluster.set$New <- NULL  # should be all FALSE
  cluster.set <- cluster.set[order(ClusterID),]
  
  # fit edge probability decay model (e.g., fit.decay="colyear")
  if (!is.na(time.var)) {
    if (!is.element(time.var, names(obj$seq.info)))
      stop(time.var, "is not a variable in obj$seq.info!")
    # fit binomial model to bipartite graph
    fit <- fit.decay(obj, time.var)
    
    # the predicted edge probability for a cluster is a function of mean node age
    #cluster.set[, "EdgeProb" := ]
  }
  
  # Attach growth info and set ID
  growth <- table(obj$seq.info[(New) & (Cluster %in% cluster.set$ClusterID), 
                               Cluster])
  cluster.set[, "Growth" := 0]
  cluster.set[ClusterID %in% as.numeric(names(growth)), 
              Growth := as.numeric(growth)]
  
  # useful when concatenating cluster sets induced by different thresholds
  cluster.set[, "DistThresh" := dist.thresh]
  cluster.set[, "SetID" := setID]
  
  class(cluster.set) <- c("cluster.set", class(cluster.set))
  
  return(cluster.set)
}


#' Fit binomial regression model to distribution of bipartite edges between
#' samples at different time points, as a model of decay in edge density with
#' time.  (internal)
#' @param edges:  data.frame, edges filtered by distance threshold, with nodes
#'                identified by integer index to time vector
#' @param times:  numeric, vector of time values for nodes
#' @return glm object
fit.decay <- function(edges, times) {
  time.counts <- table(times)
  if (sum(time.counts == 1) / length(time.counts) > 0.5) {
    stop("fit.decay only supports discrete time, ", time.var, 
         " appears to be continuous")
  }
  
  # remove new cases and associated edges
  max.time <- max(times)
  keep <- which(times < max.time)
  old.edges <- edges[edges$ID1 %in% keep & edges$ID2 %in% keep, ]
  old.edges$t1 <- times[old.edges$ID1]
  old.edges$t2 <- times[old.edges$ID2]
  
  # for every node, find the shortest edge from an older node
  # FIXME: this could be done once only outside this function
  positives <- lapply(keep, function(child) {
    child.t <- times[child]
    my.edges <- old.edges[(old.edges$ID1==child | old.edges$ID2==child), ]
    parents <- ifelse(my.edges$ID1==child, my.edges$ID2, my.edges$ID1)
    parents.t <- times[parents]
    
    my.edges$dt <- child.t - parents.t
    retro.edges <- my.edges[parents.t < child.t, ]
    retro.edges$dt[which.min(retro.edges$Distance)]
  })
  positives <- unlist(positives)
  
  # prepare data frame
  time.vals <- sort(unique(times[keep]), decreasing=TRUE)
  dts <- apply(expand.grid(time.vals, time.vals), 1, diff)
  dts <- unique(dts[dts>0])
  counts <- data.frame(dt=dts)
  counts$positives <- as.integer(table(factor(positives, levels=dts)))
  
  # calculate denominators
  counts$total <- 0
  for (i in 1:(length(time.vals)-1)) {
    t2 <- time.vals[i]
    for (j in (i+1):length(time.vals)) {
      t1 <- time.vals[j]
      dt <- t2 - t1
      possible.edges <- as.integer(
        time.counts[as.character(t1)] * time.counts[as.character(t2)]
      )
      counts$total[counts$dt==dt] <- counts$total[counts$dt==dt] + 
        possible.edges
    }
  }
  
  fit <- glm(cbind(positives, total) ~ as.integer(dt), family="binomial", 
             data=counts)
  return(fit)
}
