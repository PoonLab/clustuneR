#' Obtain paraphyletic clusters
#'
#' Clusters are defined as a group of tips diverging from a high confidence common 
#' ancestor. Divergence must be done through a series of short branches, which 
#' the branch.thresh constrains.
#'
#' @param obj: S3 object of class "phylo". The input tree file, extended to be annotated with vertex, edge and 
#' growth information.
#' @param branch.thresh: The maximum branch length criterion defining clusters.
#' A branch exceeding this value separates the tip from it's ancestor's cluster.
#' Higher values imply larger average cluster sizes.
#' @param boot.thresh: The minimum bootstrap criterion defining clusters.
#' Lower values imply larger average cluster size
#' @param setID: A numeric identifier for this cluster set.
#' @return: A set of clusters as a data.table. See example cluster.ex object 
#' documentation for an example of clustered sequence data + meta data
#' @export
#' @example examples/step.cluster_ex.R
step.cluster <- function(obj, branch.thresh = 0.03, boot.thresh = 0, setID = 0) {

  # Input Checking
  if (!is.numeric(branch.thresh) | !is.numeric(boot.thresh)) {
    stop("Clustering criteria must be numeric values")
  }
  if (!("path.info" %in% names(obj))) {
    stop("path.info must be defined for input tree `obj`, did you run ",
         "extend.tree()?")
  }
  if (!("growth.info" %in% names(obj))) {
    stop("growth.info must be defined for input tree `obj`, did you run ",
         "extend.tree()?")
  }

  # For each node (i) in the tree, find the branch on the path from (i) to 
  # the root at which the total path length exceeds the threshold.
  path.stop <- lapply(obj$node.info$Paths, function(p) {
    bl <- obj$node.info$BranchLength[p]
    cml.bl <- cumsum(bl)
    ht <- which(cml.bl > branch.thresh)[1]  # height, index into other vectors
    boots <- obj$node.info$Bootstrap[p]
    return(c(Node=p[ht], Boot=boots[ht], BranchLength=bl[ht], Height=ht))
  })
  path.stop <- as.data.frame(do.call(rbind, path.stop))
  # handle any paths that hit root before threshold length
  path.stop$Node[is.na(path.stop$Node)] <- Ntip(obj) + 1  # root index
  
  # Check bootstrap requirements, stepping back down clustered paths until 
  # they're met.
  low.support <- which(path.stop$Boot < boot.thresh)
  for (i in low.support) {
    ht <- path.stop$Height[i]
    this.path <- obj$node.info$Paths[[i]]
    boots <- obj$node.info$Bootstrap[this.path]
    blens <- obj$node.info$BranchLength[this.path]
    new.idx <- which(boots[ht:length(boots)] >= boot.thresh)[1]
    new.ht <- ht+new.idx-1
    if (is.na(new.idx)) {
      # FIXME: this should never happen because root bootstrap is set to max
      stop("step.clustering: no ancestral nodes met bootstrap threshold.")
    }
    path.stop$Node[i] <- this.path[new.ht]
    path.stop$Boot[i] <- boots[new.ht]
    path.stop$BranchLength[i] <- blens[new.ht]
    path.stop$Height[i] <- new.ht
  }
  

  # assign cluster memberships for all nodes in tree (not include new tips)
  seq.cols <- colnames(obj$seq.info)
  obj$node.info$Cluster <- path.stop["Node", ]
  
  # cluster assignments for tips only (including "new" sequences)
  obj$seq.info$Cluster <- 0
  obj$seq.info$Cluster[!obj$seq.info$New] <- obj$node.info$Cluster[1:Ntip(obj)]
  
  # build a data table of known cases (i.e., not new cases)
  cluster.set <- obj$seq.info[!(New), lapply(seq.cols, function(nm) {
    list(get(nm))
  }), by = Cluster]
  cluster.set <- cluster.set[order(Cluster), ]
  
  # collect descendants for each known case
  des <- obj$node.info[
    Cluster %in% cluster.set$Cluster, list(.(NodeID)), by = Cluster
    ]
  des <- des[order(Cluster), ]
  cluster.set[, "Descendants" := des$V1]
  cluster.set[, "Size" := length(V1[[1]]), by = 1:nrow(cluster.set)]
  colnames(cluster.set) <- c("ClusterID", seq.cols, "Descendants", "Size")
  cluster.set$New <- NULL
  
  # Assign growth cases to clusters, summing certainty for each
  obj$growth.info[, "Cluster" := obj$node.info[
    obj$growth.info$NeighbourNode, Cluster]
    ]
  obj$growth.info[(TermDistance) >= branch.thresh, Cluster := NA]
  
  growth <- obj$growth.info[!is.na(Cluster), sum(Bootstrap), by = .(Header, Cluster)]
  if (length(growth)>0){
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


#' Obtain monophyletic clusters based on pairwise distances
#'
#' Clusters are defined as  as a monophyletic clade under a high-confidence common 
#' ancestor. Some measure of divergence (criterion) within this clade must fall 
#' under a distance threshold in order for it to be labelled a cluster.
#'
#' @param t: The input tree file, annotated with vertex and edge information
#' @param dist.criterion: A particular column in node.info that must be less than 
#' a distance threshold. By default, this is "max.patristic.dist", however, the 
#' "mean.patristic.dist" is also a column. Other columns added to node.info can 
#'  be checked against dist.thresh, however, these would be added by the user.
#' @param dist.thresh: The threshold required for clustering. Monophyletic groups
#' with criterion under this value could be considered clusters.
#' @param setID: A numeric identifier for this cluster set.
#' @param boot.thresh: The minimum bootstrap threshold defining clusters. Monophyletic
#' groups with a parent node more certain than this criterion could be considered 
#' clusters
#' @return: A set of clusters as a data.table. See example cluster.ex object 
#' documentation for an example of clustered sequence data + meta data
mono.pat.cluster <- function(t, dist.thresh, boot.thresh = 0, 
                             dist.criterion = "max.patristic.dist", setID = 0) {

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


