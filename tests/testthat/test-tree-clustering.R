test_that("assign.sstrees works", {
  # prepare test fixture
  seqs <- ape::read.FASTA(test_path("test.fasta"))
  seq.info <- parse.headers(names(seqs), 
                            var.names=c("accn", "coldate", "subtype"),
                            var.transformations = c(as.character, as.Date, as.factor))
  phy <- ape::read.tree(test_path("test-oldseq.fasta.treefile"))
  phy <- import.tree(phy, seq.info = seq.info)
  log.file <- test_path("test-oldseq.fasta.log")
  phy <- extend.tree(phy, seqs, log.file)
  
  # bypass bootstrap threshold
  result <- assign.sstrees(phy, branch.thresh=0.04, boot.thresh=0)
  expect_true(is.data.frame(result))
  expect_equal(names(result), c("Node", "Boot", "BranchLength", "Height"))
  
  expected <- c(8, 8, 10, 4, 11, 11, 7, 9, 7, 7, 11)
  expect_equal(result$Node, expected)
  
  expected <- c(0.99, 0.99, 0.90, 1, 1, 1, 1, 1, 1, 1, 1)
  expect_equal(result$Boot, expected)
  
  expected <- c(0.036096058, 0.036096058, 0.025737776, 0.074902943, 
                0.072356702, 0.072356702, NA, 0.004451792, NA, NA, 0.072356702)
  expect_equal(result$BranchLength, expected)
  
  expected <- c(2, 2, 2, 1, 2, 2, 1, 2, 2, 3, 1)
  expect_equal(result$Height, expected)
  
  # reducing the path length threshold should change the cluster assignments
  result <- assign.sstrees(phy, branch.thresh=0.03, boot.thresh=0)
  expected <- c(8, 2, 3, 4, 5, 11, 7, 8, 7, 9, 11)
  expect_equal(result$Node, expected)
  
  # now check bootstrap support
  result <- assign.sstrees(phy, branch.thresh=0.03, boot.thresh=1)
  expected <- c(9, 2, 3, 4, 5, 11, 7, 9, 7, 9, 11)
  expect_equal(result$Node, expected)
  
  result <- assign.sstrees(phy, branch.thresh=0.03, boot.thresh=0.98)
  expected <- c(8, 2, 3, 4, 5, 11, 7, 8, 7, 9, 11)
  expect_equal(result$Node, expected)
})

test_that("step.cluster works", {
  # prepare test fixture
  seqs <- ape::read.FASTA(test_path("test2.fasta"))
  seq.info <- parse.headers(names(seqs), 
                            var.names=c("accn", "coldate", "subtype"),
                            var.transformations = c(as.character, as.Date, as.factor))
  phy <- ape::read.tree(test_path("test2-old.fasta.treefile"))
  phy <- import.tree(phy, seq.info = seq.info)
  log.file <- test_path("test2-old.fasta.log")
  phy <- extend.tree(phy, seqs, log.file)
  
  result <- step.cluster(phy, branch.thresh=0.04, boot.thresh=0)
  
})
