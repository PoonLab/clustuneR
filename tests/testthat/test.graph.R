source("R/graph.setup.R")
source("R/graph.clustering.R")
source("R/analysis.R")
library(testthat)
library(ape)

load("data/seq-phylo_ex.RData")
load("data/seq.info_ex.RData")
load("data/graphs_ex.RData")

edge.info.tn93 <- dist.dna(seqs.full, pairwise.deletion = T, as.matrix = T, model = "TN93", )
edge.info.patristic <- cophenetic.phylo(tree.full)

test_that("Graph can be created from dist.dna TN93 with no new seqs", {
  expect_warning(expect_error(create.graph(seq.info, edge.info.tn93), NA),
                 "No new sequences are specified by a New column in seq.info.")
})

test_that("Graph can be created from cophenetic.phylo pairwise patristic with no new seqs", {
  expect_warning(expect_error(create.graph(seq.info, edge.info.patristic), NA),
                 "No new sequences are specified by a New column in seq.info.")
})

test_that("Component clusters can be defined at boundary thresholds with no growth", {
  param.list <- lapply(c(-Inf,0.007, Inf), function(x){list("g"=g.nogrowth, "dist.thresh"=x)})
  expect_error(multi.cluster(component.cluster, param.list, mc.cores = 4), NA)
})

test_that("Component clusters can be defined at boundary thresholds with growth", {
  param.list <- lapply(c(-Inf,0.007, Inf), function(x){list("g"=g.growing, "dist.thresh"=x)})
  expect_error(multi.cluster(component.cluster, param.list, mc.cores = 4), NA)
})
