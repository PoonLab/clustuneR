setwd("../../")

source("R/tree.setup.R")
source("R/sequence.setup.R")
source("R/tree.clustering.R")
source("R/analysis.R")
library(testthat)
library(ape)
library(phangorn)

load("data/seq.info_ex.RData")
load("data/seq-phylo_ex.RData")
load("data/tree_ex.RData")

tree.old.sample <- drop.tip(tree.old, 1:30)
seqs.old <- seqs.full[which(names(seqs.full)%in%tree.old.sample$tip.label)]
seqs.new <- seqs.full[which(!(names(seqs.full)%in%tree.old.sample$tip.label))]
seqs.new.sample <- seqs.new[which(names(seqs.new)%in%seq.info$Header)[1:3]]
full.align <- c(seqs.old, seqs.new.sample)

test_that("Tree can be extended with no new seqs", {
  expect_warning(expect_error(extend.tree(t=tree.old.sample, seq.info=seq.info, mc.cores=1), NA),
                 "Ignoring growth information, path to logfile and full sequence alignment required.")
})

test_that("IQ-TREE  Logfiles can be recognized and translated", {
  expect_silent(translate.log("data/IQTREE_log_ex.txt"))
})

test_that("RAxML Logfiles can be recognized and translated", {
  expect_silent(translate.log("data/RAxML_log_ex.txt"))
})

param.list <- lapply(c(-Inf,0.03, Inf), function(x){list("t"=t, "branch.thresh"=x)})
clusters <- multi.cluster(step.cluster, param.list, mc.cores = 1)

test_that("Total growth increases with more liberal thresholds", {
  expect_equal(max(clusters[BranchThresh==-Inf, Growth]), 0)
  expect_equal(min(clusters[BranchThresh==Inf, Growth]), nrow(t$seq.info[(New),]))
  expect_equal(max(clusters[BranchThresh==Inf, Growth]), max(clusters$Growth))
})

test_that("Clusters range from completely separated to completely agglomerated", {
  expect_equal(nrow(clusters[BranchThresh==-Inf]), nrow(t$seq.info[!(New),]))
  expect_equal(nrow(clusters[BranchThresh==Inf]), 1)
})
