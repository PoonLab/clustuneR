source("R/tree.setup.R")
source("R/sequence.setup.R")
source("R/tree.clustering.R")
source("R/analysis.R")
library(testthat)
library(ape)
library(phangorn)

load("data/seq.info_ex.RData")
load("data/seq-phylo_ex.RData")
load("data/trees_ex.RData")

test_that("Tree can be extended with no new seqs", {
  expect_warning(expect_error(extend.tree(t=tree.full, seq.info=seq.info, mc.cores=4), NA),
                 "Ignoring growth information, path to logfile and full sequence alignment required.")
})

test_that("Tree can be extended with new seqs to represent growth", {
  expect_warning(expect_error(extend.tree(tree.old, seq.info, mc.cores=4, full.align = seqs.full, 
                              log.file = "data/IQTREE_log_ex.txt"), NA),
                 "Not all newly added sequences are noted in the seq.info of the tree")
})

test_that("FastTree Logfiles can be recognized and translated", {
  expect_silent(translate.log("data/FastTree_log_ex.txt"))
})

test_that("IQ-TREE  Logfiles can be recognized and translated", {
  expect_silent(translate.log("data/IQTREE_log_ex.txt"))
})

test_that("RAxML Logfiles can be recognized and translated", {
  expect_silent(translate.log("data/RAxML_log_ex.txt"))
})

test_that("Step clusters can be defined at boundary thresholds with no growth", {
  param.list <- lapply(c(-Inf,0.007, Inf), function(x){list("t"=t.nogrowth, "branch.thresh"=x)})
  expect_error(multi.cluster(step.cluster, param.list, mc.cores = 4), NA)
})

test_that("Step clusters can be defined at boundary thresholds with growth", {
  param.list <- lapply(c(-Inf,0.007, Inf), function(x){list("t"=t.growing, "branch.thresh"=x)})
  expect_error(multi.cluster(step.cluster, param.list, mc.cores = 4), NA)
})

test_that("Monophyletic patristic clusters can be defined at boundary thresholds with no growth", {
  param.list <- lapply(c(-Inf,0.007, Inf), function(x){list("t"=t.nogrowth, "dist.thresh"=x)})
  expect_error(multi.cluster(mono.pat.cluster, param.list, mc.cores = 4), NA)
})