test_that("Graph can be created from dist.dna and cophenetic.phylo TN93 with no new seqs", {
  edge.info.tn93 <- ape::dist.dna(alignment.ex, pairwise.deletion = T, as.matrix = T, model = "TN93", )
  edge.info.patristic <- ape::cophenetic.phylo(full.tree.ex)

  expect_error(create.graph(seq.info.ex, edge.info.tn93), NA)
  expect_error(create.graph(seq.info.ex, edge.info.patristic), NA)
})

param.list <- lapply(c(-Inf,0.06, Inf), function(x){list("g"=graph.ex, "dist.thresh"=x)})
clusters <- multi.cluster(component.cluster, param.list, mc.cores = 1)

test_that("Total growth increases with more liberal thresholds", {
  expect_equal(max(clusters[DistThresh==-Inf, Growth]), 0)
  expect_equal(min(clusters[DistThresh==Inf, Growth]), nrow(graph.ex$seq.info[(New),]))
  expect_equal(max(clusters[DistThresh==Inf, Growth]), max(clusters$Growth))
})

test_that("Clusters range from completely separated to completely agglomerated", {
  expect_equal(nrow(clusters[DistThresh==-Inf]), nrow(graph.ex$seq.info[!(New),]))
  expect_equal(nrow(clusters[DistThresh==Inf]), 1)
})

