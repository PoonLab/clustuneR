#'Generate data found in /data folder
#'
#'This is partially intended as example use code, however may also act as secondary,
#'informal testing in the development cycle and as a tool to update data quickly if required.
generate.all <- function() {

  generate.seq.info()
  generate.graphs()
  generate.trees()
  generate.clusters()
  generate.fit.results()

}

#' Obtain basic sequence information
generate.seq.info <- function() {

  load("data/seq-phylo_ex.RData")

  seq.info <- pull.headers(seqs.full,var.names = c("ID", "CollectionDate", "Subtype"),
                           var.transformations =list(as.character, as.Date, as.factor))

  save("seq.info", file="data/seq.info_ex.RData")

}

#' Create example graphs
generate.graphs <- function() {

  load("data/seq-phylo_ex.RData")
  load("data/seq.info_ex.RData")

  edge.info <- ape::dist.dna(seqs.full, pairwise.deletion = T, as.matrix = T, model = "TN93", )
  g.nogrowth <- create.graph(seq.info, edge.info)

  new.year <- max(seq.info$CollectionDate) - 365
  which.new <- which(seq.info$CollectionDate > new.year)
  annotate.new(seq.info,which.new)
  g.growing <- create.graph(seq.info, edge.info)

  save(list=c("g.nogrowth","g.growing"), file="data/graphs_ex.RData")
}

#' Create example trees
generate.trees <- function(){

  load("data/seq-phylo_ex.RData")
  load("data/seq.info_ex.RData")

  t.nogrowth <- extend.tree(tree.full, seq.info, mc.cores=4)
  t.growing <- extend.tree(tree.old, seq.info, mc.cores=4,
                           full.align = seqs.full, log.file = "data/IQTREE_log_ex.txt")
  save(list=c("t.nogrowth","t.growing"), file="data/trees_ex.RData")
}

#' Create example clusters
generate.clusters <- function() {

  load("data/graphs_ex.RData")
  load("data/trees_ex.RData")

  param.list <- lapply(c(-Inf,0.007, Inf), function(x){list("t"=t.growing, "branch.thresh"=x)})
  step.cluster.data <- multi.cluster(step.cluster, param.list, mc.cores = 4)

  param.list <- lapply(c(-Inf,0.007, Inf), function(x){list("t"=t.nogrowth, "dist.thresh"=x)})
  mono.pat.cluster.data <- multi.cluster(mono.pat.cluster, param.list, mc.cores = 4)

  param.list <- lapply(c(-Inf,0.007, Inf), function(x){list("g"=g.growing, "dist.thresh"=x)})
  component.cluster.data <- multi.cluster(component.cluster, param.list, mc.cores = 4)

  component.cluster.set <- component.cluster(g.growing, dist.thresh = 0.015)

  save(list=c("component.cluster.data","step.cluster.data", "mono.pat.cluster.data", "component.cluster.set"), file="data/clusters_ex.RData")

}

#' Create analysis from a set of clusters
generate.fit.results <- function() {

  load("data/clusters_ex.RData")

  predictive.models <- list(
    "NullModel" = function(x){
      glm(Growth~Size, data=x, family="poisson")
    },
    "TimeModel" = function(x){
      glm(Growth~Size+CollectionDate, data=x, family="poisson")
    }
  )

  predictor.transformations <- list(
    "CollectionDate" = function(x){mean(x)}
  )

  param.list.boot <- lapply(seq(0,0.04,0.001), function(x){list("t"=t.growing, "branch.thresh"=x, "boot.thresh"=0.9)})
  param.list.noboot <- lapply(seq(0,0.04,0.001), function(x){list("t"=t.growing, "branch.thresh"=x)})

  step.cluster.data.boot <- multi.cluster(step.cluster, param.list.boot, mc.cores = 4)
  step.cluster.data.noboot <- multi.cluster(step.cluster, param.list.noboot, mc.cores = 4)

  res.step.boot <- fit.analysis(step.cluster.data.boot,
                                predictor.transformations = predictor.transformations,
                                predictive.models=predictive.models)
  res.step.noboot <- fit.analysis(step.cluster.data.noboot,
                                  predictor.transformations = predictor.transformations,
                                  predictive.models=predictive.models)

  save(list=c("res.step.boot", "res.step.noboot"), file="data/fit.results_ex.RData")
}
