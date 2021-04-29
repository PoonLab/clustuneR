#'Generate data found in /data folder
#'
#'This is partially intended as example use code, however may also act as informal 
#'testing in the development cycle and as a tool to update data quickly if required.
generate.all <- function() {
  generate.seq.info()
  generate.graph()
  generate.tree()
}

#' Obtain basic sequence information
generate.seq.info <- function() {

  seq.info.ex <- pull.headers(alignment.ex, var.names = c("ID", "CollectionDate", "Subtype"),
                           var.transformations =list(as.character, as.Date, as.factor))

  save(seq.info.ex, file="data/seq.info.ex.RData")

}

#' Create example graphs
generate.graph <- function() {

  edge.info <- ape::dist.dna(alignment.ex, pairwise.deletion = T, as.matrix = T, model = "TN93", )

  new.year <- max(seq.info.ex$CollectionDate) - 365
  which.new <- which(seq.info.ex$CollectionDate > new.year)
  graph.ex <- create.graph(seq.info.ex, edge.info, which.new)

  save(graph.ex, file="data/graph.ex.RData")
}

#' Create example trees
generate.tree <- function(){

  extended.tree.ex <- extend.tree(old.tree.ex, seq.info.ex, full.align = alignment.ex, log.file = "data/IQTREE_log_ex.txt")

  save(extended.tree.ex, file="data/extended.tree.ex.RData")
}

#' Create example cluster set
generate.cluster <- function(){

  cluster.ex <- component.cluster.set <- component.cluster(graph.ex, dist.thresh = 0.055)

  save(cluster.ex, file="data/cluster.ex.RData")
}
