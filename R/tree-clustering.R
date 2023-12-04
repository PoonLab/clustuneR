#' Obtain paraphyletic clusters
#'
#' Clusters are defined as a group of tips diverging from a high confidence 
#' common ancestor. Divergence must be done through a series of short branches, 
#' which the branch.thresh constrains.
#'
#' @param obj: S3 object of class "phylo". The input tree file, extended to be 
#'             annotated with vertex, edge and growth information.
#' @param branch.thresh: numeric, the maximum branch length criterion defining 
#'                       clusters.  A branch exceeding this value separates the 
#'                       tip from its ancestor's cluster.  Higher values imply 
#'                       larger average cluster sizes.
#' @param boot.thresh: numeric, the minimum bootstrap criterion defining 
#'                     clusters.  Lower values imply larger average cluster 
#'                     size.
#' @param setID:  integer, a numeric identifier for this cluster set.
#' @param resolve:  character, if 'max', assign new case to cluster with highest
#'                  probability; if 'random', assign at random in proportion 
#'                  to probability of placement.
#' @return:  A set of clusters as a data.table.
#' @export
step.cluster <- function(obj, branch.thresh=0.03, boot.thresh=0, setID=0, 
                         resolve='max') {
  if (!is.numeric(branch.thresh) | !is.numeric(boot.thresh)) {
    stop("Clustering criteria must be numeric values")
  }
  if (!("growth.info" %in% names(obj))) {
    stop("growth.info must be defined for input tree `phy`, did you run ",
         "extend.tree()?")
  }

  # assign cluster memberships for all nodes in tree (not include new tips)
  phy <- assign.sstrees(obj, branch.thresh, boot.thresh)
  
  # build a data table of known cases (i.e., not new cases)
  seq.cols <- colnames(phy$seq.info)
  cluster.set <- phy$seq.info[!(New), lapply(seq.cols, function(nm) {
    list(get(nm))
  }), by = Cluster]
  cluster.set <- cluster.set[order(Cluster), ]
  names(cluster.set) <- c("Cluster", seq.cols)
  cluster.set$New <- NULL
  
  # collect descendants for each known case to calculate cluster sizes
  des <- sapply(
    split(phy$node.info$Descendants[1:Ntip(phy)], 
          phy$node.info$Cluster[1:Ntip(phy)]), 
    function(x) unique(unlist(x))
    )
  cluster.set[, "Descendants" := des]
  cluster.set$Size <- sapply(cluster.set$Descendants, length)
  
  # Assign growth cases to clusters, sum over branch placements within clusters
  growth <- phy$growth.info[!is.na(Cluster), sum(Bootstrap), 
                            by=.(Header, Cluster)]
  if (resolve=='max') {
    growth <- growth[, Cluster[which.max(V1)], by=.(Header)]
  } 
  else if (resolve=='random') {
    growth <- growth[, Cluster[sample(1:length(V1), size=1, prob=V1)], 
                     by=.(Header)]
  }
  else {
    stop("Unrecognized `resolve` type", resolve, ". Expected `max` or `random`.")
  }
  
  growth <- table(growth$V1)
  # FIXME: is this really necessary?
  growth <- growth[which(as.numeric(names(growth)) %in% cluster.set$Cluster)]
  
  # Attach growth info and a set ID to clusters
  cluster.set[, "Growth" := 0]
  cluster.set[Cluster %in% as.numeric(names(growth)), 
              "Growth":=as.numeric(growth)]
  
  cluster.set[, "BranchThresh" := branch.thresh]
  cluster.set[, "BootThresh" := boot.thresh]
  cluster.set[, "SetID" := setID]

  return(cluster.set)
}


#' Assign subset trees
#' 
#' For each node (i) in the tree, find the branch on the path from (i) to 
#' the ancestral node at which the total path length exceeds the threshold.
#' If bootstrap support at that node is below the threshold, step further
#' down the tree (toward global root) until  it is met.
#' @param phy: an object of class ape::phylo.
#' @param branch.thresh:  numeric, branch length threshold (distance from
#'                        node to root of subset tree).
#' @param boot.thresh:  numeric, bootstrap support threshold.
#' @param debug: logical, use TRUE for unit testing (expose data frame)
#' @return data frame with a row for each node in the tree, identifying 
#'         the subset tree-defining root node
#' @export
assign.sstrees <- function(phy, branch.thresh, boot.thresh, debug=FALSE) {
  res <- lapply(phy$node.info$Paths, function(p) {
    boots <- phy$node.info$Bootstrap[p]
    bl <- phy$node.info$BranchLength[p]
    cml.bl <- cumsum(bl)
    
    ht <- which(cml.bl > branch.thresh)[1]  # height, index into other vectors
    if (is.na(ht)) ht <- length(p)  # past root of tree
    return(c(Node=p[ht], Boot=boots[ht], BranchLength=bl[ht], Height=ht))
  })
  res <- as.data.frame(do.call(rbind, res))
  
  # Check bootstrap requirements, stepping back down clustered paths until 
  # they're met.
  low.support <- which(res$Boot < boot.thresh)
  for (i in low.support) {
    ht <- res$Height[i]  # current position in path for this node
    
    # extract attributes for all nodes on this path
    this.path <- phy$node.info$Paths[[i]]
    boots <- phy$node.info$Bootstrap[this.path]
    blens <- phy$node.info$BranchLength[this.path]

    # determine new stopping point, bypassing path length threshold
    new.idx <- which(boots[ht:length(boots)] >= boot.thresh)[1]
    new.ht <- ht+new.idx-1
    if (is.na(new.idx)) {
      # this should never happen because root bootstrap is set to max
      stop("assign.sstrees: no ancestral nodes met bootstrap threshold.")
    }
    
    # update data frame with new entry
    res$Node[i] <- this.path[new.ht]
    res$Boot[i] <- boots[new.ht]
    res$BranchLength[i] <- blens[new.ht]
    res$Height[i] <- new.ht
  }
  
  if(debug) return(res)  # for unit testing
  
  # transfer cluster assignments to phylo object
  phy$node.info$Cluster <- res$Node
  
  # cluster assignments for tips only (including "new" sequences)
  phy$seq.info$Cluster <- 0
  phy$seq.info$Cluster[!phy$seq.info$New] <- phy$node.info$Cluster[1:Ntip(phy)]
 
  phy$growth.info[, "Cluster" := phy$node.info[
    phy$growth.info$NeighbourNode, Cluster]
  ]
  # new cases that are too far from any known case are not in clusters
  phy$growth.info[(TermDistance) >= branch.thresh, Cluster := NA]
  
  return(phy)
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


