test_that("Tree can be extended with no new seqs, throws warning", {
  expect_warning(expect_error(extend.tree(t=old.tree.ex, seq.info=seq.info.ex, mc.cores=1), NA),
                 "Ignoring growth information, path to logfile and full sequence alignment required.")
})

test_that("IQ-TREE  Logfiles can be recognized and translated", {
  expect_silent(translate.log("../../data/IQTREE_log_ex.txt"))
})

test_that("RAxML Logfiles can be recognized and translated", {
  expect_silent(translate.log("../../data/RAxML_log_ex.txt"))
})

param.list <- lapply(c(-Inf,0.03, Inf), function(x){list("t"=extended.tree.ex, "branch.thresh"=x)})
clusters <- multi.cluster(step.cluster, param.list, mc.cores = 1)

test_that("Total growth increases with more liberal thresholds", {
  expect_equal(max(clusters[BranchThresh==-Inf, Growth]), 0)
  expect_equal(min(clusters[BranchThresh==Inf, Growth]), nrow(extended.tree.ex$seq.info[(New),]))
  expect_equal(max(clusters[BranchThresh==Inf, Growth]), max(clusters$Growth))
})

test_that("Clusters range from completely separated to completely agglomerated", {
  expect_equal(nrow(clusters[BranchThresh==-Inf,]), nrow(extended.tree.ex$seq.info[!(New),]))
  expect_equal(nrow(clusters[BranchThresh==Inf,]), 1)
})
