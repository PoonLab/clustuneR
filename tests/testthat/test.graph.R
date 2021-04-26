setwd("../../")

source("R/graph.setup.R")
source("R/graph.clustering.R")
source("R/analysis.R")
library(testthat)
library(ape)

load("data/seq-phylo_ex.RData")
load("data/seq.info_ex.RData")
load("data/graph_ex.RData")



test_that("Graph can be created from dist.dna and cophenetic.phylo TN93 with no new seqs", {
  edge.info.tn93 <- dist.dna(seqs.full, pairwise.deletion = T, as.matrix = T, model = "TN93", )
  edge.info.patristic <- cophenetic.phylo(tree.full)

  expect_error(create.graph(seq.info, edge.info.tn93), NA)
  expect_error(create.graph(seq.info, edge.info.patristic), NA)
})

param.list <- lapply(c(-Inf,0.06, Inf), function(x){list("g"=g, "dist.thresh"=x)})
clusters <- multi.cluster(component.cluster, param.list, mc.cores = 1)

test_that("Total growth increases with more liberal thresholds", {
  expect_equal(max(clusters[DistThresh==-Inf, Growth]), 0)
  expect_equal(min(clusters[DistThresh==Inf, Growth]), nrow(g$seq.info[(New),]))
  expect_equal(max(clusters[DistThresh==Inf, Growth]), max(clusters$Growth))
})

test_that("Clusters range from completely separated to completely agglomerated", {
  expect_equal(nrow(clusters[DistThresh==-Inf]), nrow(g$seq.info[!(New),]))
  expect_equal(nrow(clusters[DistThresh==Inf]), 1)
})

