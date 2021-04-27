#' Obtain clusters based on path to ancestor in tree
#'
#' Clusters are defined as a series of tips diverging from A high confidence common ancestor.
#' This divergence must be done through a series of short branches, which the branch.thresh constrains.
#'
#' @param t: The input tree file, annotated with vertex and edge information
#' @param branch.thresh: The maximum branch length criterion defining clusters
#' @param boot.thresh: The minimum bootstrap criterion defining clusters
#' @param setID: If several different parameter ranges are used, the setID can identify them
#' @return A data table which extends a subset of node.info. This includes growth info
#' @export
#' @example examples/step.cluster_ex.R
step.cluster <- function(t, branch.thresh = 0.03, boot.thresh = 0, setID = 0) {

  # Input Checking
  if (!is.numeric(branch.thresh) | !is.numeric(boot.thresh)) {
    stop("Clustering criteria must be numeric values")
  }
  if (!("path.info" %in% names(t))) {
    stop("path.info must be defined for tree")
  }
  if (!("growth.info" %in% names(t))) {
    stop("growth.info must be defined for tree")
  }

  # Obtain the stopping point in the path based on branch.thresh
  path.stop <- sapply(t$path.info, function(p) {
    h <- which(p["BranchLength", ] > branch.thresh)[1]
    c(p[, h], h)
  })
  rownames(path.stop)[4] <- "Height"
  path.stop["Node", is.na(path.stop["Node", ])] <- length(t$tip.label) + 1

  # Check bootstrap requirements, stepping back down clustered paths until they're met.
  i <- which(path.stop["Boot", ] < boot.thresh)
  if (length(i) > 0) {
    path.stop[, i] <- sapply(i, function(j) {
      p <- t$path.info[[j]]
      p.boots <- p["Boot", 1:path.stop["Height", j]]
      new.h <- which(p.boots >= boot.thresh)[1]
      return(c(p[, new.h], new.h))
    })
  }

  # Assign Clusters and update membership info for each
  seq.cols <- colnames(t$seq.info)
  t$node.info[, "Cluster"] <- path.stop["Node", ]
  t$seq.info[, "Cluster" := 0]
  t$seq.info[!(New), Cluster := t$node.info[1:length(t$tip.label), Cluster]]

  cluster.set <- t$seq.info[!(New), lapply(seq.cols, function(nm) {
    list(get(nm))
  }), by = Cluster]
  cluster.set <- cluster.set[order(Cluster), ]
  des <- t$node.info[Cluster %in% cluster.set$Cluster, list(.(NodeID)), by = Cluster]
  des <- des[order(Cluster), ]
  cluster.set[, "Descendants" := des$V1]
  cluster.set[, "Size" := length(V1[[1]]), by = 1:nrow(cluster.set)]
  colnames(cluster.set) <- c("ClusterID", seq.cols, "Descendants", "Size")
  cluster.set$New <- NULL

  # Assign growth cases to clusters, summing certainty for each
  t$growth.info[, "Cluster" := t$node.info[t$growth.info$NeighbourNode, Cluster]]
  t$growth.info[(TermDistance) >= branch.thresh, Cluster := NA]

  growth <- t$growth.info[!is.na(Cluster), sum(Bootstrap), by = .(Header, Cluster)]
  if(length(growth)>0){
    growth <- growth[V1 >= boot.thresh, Cluster[which.max(V1)], by = .(Header)]
  }
  growth <- table(growth$V1)
  growth <- growth[which(as.numeric(names(growth)) %in% cluster.set$ClusterID)]

  # Attach growth info and a set ID to clusters
  cluster.set[, "Growth" := 0]
  cluster.set[ClusterID %in% as.numeric(names(growth)), "Growth" := as.numeric(growth)]

  cluster.set[, "BranchThresh" := branch.thresh]
  cluster.set[, "BootThresh" := boot.thresh]
  cluster.set[, "SetID" := setID]

  return(cluster.set)
}


## - TO-DO: SOLVE MONOPHYLETIC CLUSTER GROWTH IN A SIMPLE WAY -##
#' Obtain clusters based on a monophyletic group in tree
#'
#' Clusters as a monophyletic clade under a high-confidence common ancestor.
#' The pairwise patristic distances in this clade must all
#'
#' @param t: The input tree file, annotated with vertex and edge information
#' @param dist.criterion: A particular column in node.info that must be less than a distance threshold
#' @param dist.thresh: The threshold required for clustering.
#' @param setID: If several different parameter ranges are used, the setID can identify them.
#' @param boot.thresh: The minimum bootstrap criterion defining clusters
#' @return A data table which extends a subset of node.info. This includes growth info
mono.pat.cluster <- function(t, dist.thresh, boot.thresh = 0, dist.criterion = "max.patristic.dist", setID = 0) {

  warning("Method unfinished. Growth information for clusters not included")

  # Input Checking
  if (!is.numeric(dist.thresh) | !is.numeric(boot.thresh)) {
    stop("Clustering criteria must be numeric values")
  }
  if (!("growth.info" %in% names(t))) {
    stop("growth.info must be defined for tree")
  }

  # Cluster Criterion checking
  t$node.info[, "Clustered" := F]
  t$node.info[(Bootstrap >= boot.thresh) & (get(dist.criterion) <= dist.thresh), "Clustered" := T]
  clustered.des <- unlist(t$node.info[(Clustered), Descendants])

  times.clustered <- t$node.info[(Clustered), length(which(clustered.des %in% NodeID)), by = which(Clustered)]
  sub.clusters <- times.clustered[((which <= length(t$tip.label)) & (V1 > 1)) | ((which > length(t$tip.label)) & (V1 > 0)), which]
  t$node.info[sub.clusters, "Clustered" := F]

  # Find parent clusters
  t$node.info[, "Cluster" := 0]
  for (i in which(t$node.info[, Clustered])) {
    t$node.info[i, "Cluster" := i]
    t$node.info[t$node.info$Descendants[[i]], "Cluster" := i]
  }
  t$seq.info[, "Cluster" := 0]
  t$seq.info[!(New), Cluster := t$node.info[1:length(t$tip.label), Cluster]]

  seq.cols <- colnames(t$seq.info)
  cluster.set <- t$seq.info[!(New), lapply(seq.cols, function(nm) {
    list(get(nm))
  }), by = Cluster]
  cluster.set[, "Size" := length(V1[[1]]), by = 1:nrow(cluster.set)]
  colnames(cluster.set) <- c("ClusterID", seq.cols, "Size")
  cluster.set$New <- NULL

  # Attach growth info and a set ID to clusters
  cluster.set[, "Growth" := 0]

  cluster.set[, "DistThresh" := dist.thresh]
  cluster.set[, "BootThresh" := boot.thresh]
  cluster.set[, "SetID" := setID]

  return(cluster.set)
}


