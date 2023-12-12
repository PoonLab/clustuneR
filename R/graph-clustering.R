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
component.cluster <- function(obj, dist.thresh, setID=0, time.var=NA) {
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
    weights <- fit.decay(
      obj, times=obj$seq.info[[time.var]], 
      dist.thresh=dist.thresh, adjusted=adjusted)
    cluster.set$Weight <- split(weights[!obj$seq.info$New], obj$seq.info$Cluster[!obj$seq.info$New])
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
#' time.
#' 
#' @param obj:  S3 object of class objEdgeData, returned from 
#' @param times:  numeric, vector of time values for nodes (from seq.info)
#' @param adjusted:  logical, if TRUE then include mean out-edge density 
#'                   as a model term to adjust for variation in sampling rates.
#' @return  numeric, weights for every node in obj$seq.info
#' @export
fit.decay <- function(obj, times, dist.thresh, adjusted=TRUE) {
  # compute some useful quantities
  edges <- obj$edge.info
  edges$t1 <- times[edges$ID1]
  edges$t2 <- times[edges$ID2]
  edges$tMax <- pmax(edges$t1, edges$t2)
  edges$tDiff <- abs(edges$t1 - edges$t2)
  
  # e-edges are filtered by distance threshold, only one in-edge per new case
  e.edges <- edges[edges$Distance < dist.thresh, ]
  time.counts <- table(times)
  if (sum(time.counts == 1) / length(time.counts) > 0.8) {
    # too many unique values
    stop("fit.decay only supports discrete time, `times` appears to be continuous")
  }
  
  # count edges that terminate in a given year
  e.times <- sapply(as.numeric(names(time.counts)), function(ti) 
    sum(e.edges$t1==ti | e.edges$t2==ti))
  names(e.times) <- names(time.counts)
  
  # f-edges are not threshold-filtered and exclude edges to new cases 
  new.cases <- which(obj$seq.info$New)
  f.edges <- edges[!(edges$ID1 %in% new.cases) & !(edges$ID2 %in% new.cases), ]

  # exclude edges between cases from the same time point
  f.edges <- f.edges[f.edges$t1 != f.edges$t2, ]
  
  # for each node, limit one edge from a previous time point
  f.edges$child <- ifelse(f.edges$tMax==f.edges$t1, f.edges$ID1, f.edges$ID2)
  temp <- lapply(split(f.edges, f.edges$child), function(x) {
    x[which.min(x$Distance), ]
  })
  f.edges <- do.call(rbind, temp)
  
  # for every time point T and preceding time point S, get:
  age.data <- lapply(unique(f.edges$tMax), function(end.t) {
    # (1) the number of nodes at time T
    total <- time.counts[as.character(end.t)]
    
    #(2) the number of edges from S to T with distance below threshold
    my.edges <- f.edges[f.edges$tMax==end.t, ]
    biparts <- split(my.edges, my.edges$tDiff)
    positives <- sapply(biparts, function(bp) sum(bp$Distance <= dist.thresh))
    
    # (3) the average number of filtered edges from a node at time S to any 
    # subsequent time point
    tdiff <- as.numeric(names(positives))
    outedge.dens <- sapply(tdiff, function(td) {
      start.t <- end.t - td
      t.key <- as.character(start.t)
      sum(e.times[t.key] / time.counts[t.key])
    })
    data.frame(positives, total, outedge.dens, tdiff)
  })
  age.data <- do.call(rbind, age.data)
  
  if (adjusted) {
    fit <- glm(cbind(positives, total) ~ tdiff+outedge.dens, data=age.data, 
               family="binomial")  
  } else {
    fit <- glm(cbind(positives, total) ~ tdiff, data=age.data, 
               family="binomial")
  }
  
  weights <- predict(
    fit, type = "response", 
    newdata = data.frame(
      tdiff = max(times) - times,
      outedge.dens = as.numeric(e.times[as.character(times)]) /
        as.numeric(time.counts[as.character(times)])
      ))
  return(weights)
}
