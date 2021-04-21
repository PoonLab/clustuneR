#' This is a testing script to run through a variety of common function uses
#' example txt log files and sequence/tree data inputs are obtained from data
load("data/test_data.RData")
devtools::load_all()

###SEQ/TREE/MODEL SETUP TESTING
seq.info <- pull.headers(seqs.full,var.names = c("ID", "CollectionDate", "Subtype"),
                         var.transformations =list(as.character, as.Date, as.factor))
seq.info$ID <- NULL
new.year <- max(seq.info$CollectionDate) - 365
which.new <- which(seq.info$CollectionDate > new.year)
seq.info <- annotate.new(seq.info,which.new)

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

#### GRAPH TESTING
edge.info.tn93 <- ape::dist.dna(seqs.full, pairwise.deletion = T, as.matrix = T, model = "TN93", )
edge.info.patristic <- ape::cophenetic.phylo(tree.full)

g.tn93 <- create.graph(seq.info, edge.info.tn93)
g.patristic <- create.graph(seq.info, edge.info.patristic)

clusters.tn93 <- component.cluster(g.tn93, 0.007)
clusters.patristic <- component.cluster(g.patristic, 0.007)

param.list.tn93 <- lapply(seq(0,0.01,0.0001), function(x){list("g"=g.tn93, "dist.thresh"=x)})
param.list.patristic <- lapply(seq(0,0.08,0.001), function(x){list("g"=g.patristic, "dist.thresh"=x)})
cluster.range.tn93 <- multi.cluster(component.cluster, param.list.tn93, mc.cores = 4)
cluster.range.patristic <- multi.cluster(component.cluster, param.list.patristic, mc.cores = 4)

res.tn93 <- fit.analysis(cluster.data=cluster.range.tn93, predictive.models=predictive.models,
                         predictor.transformations = predictor.transformations)
res.tn93 <- cbind(res.tn93, get.AIC(res.tn93))

res.patristic <- fit.analysis(cluster.data=cluster.range.patristic, predictive.models=predictive.models,
                              predictor.transformations = predictor.transformations)
res.patristic <- cbind(res.patristic, get.AIC(res.patristic))

#### PPLACER/TREE TESTING
stats.json.ft.test <- translate.log(log.file = "data/FastTree_LogEx.txt", program = "FastTree")
stats.json.rml.test <- translate.log(log.file = "data/RAxML_LogEx.txt", program = "RAxML")
stats.json <- translate.log(log.file = "data/IQTREE_old.tree_Log.txt", program = "IQ-TREE")
refpkg <- taxit.create(tree.old, seqs.full, stats.json)
tree.grown  <- run.pplacer_guppy(refpkg)
mc.cores <- 4

tree.extended <- extend.tree(tree.old, seq.info, mc.cores=mc.cores)
tree.extended$growth.info <- annotate.growth(tree.extended, tree.grown, mc.cores=mc.cores)

clusters.step <- step.cluster(tree.extended, 0.007, 0.3)
clusters.mono <- mono.pat.cluster(tree.extended, 0.07, 0.3)

param.list <- lapply(seq(0,0.04,0.001), function(x){list("t"=tree.extended, "branch.thresh"=x)})
cluster.data <- multi.cluster(step.cluster, param.list, mc.cores = 4)


res <- fit.analysis(cluster.data,
                    predictor.transformations = predictor.transformations,
                    predictive.models=predictive.models)
res <- cbind(res, get.AIC(res))





