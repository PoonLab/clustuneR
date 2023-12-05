test_that("component.cluster works", {
  # load test fixture
  edge.info <- read.csv(test_path("test.tn93.csv"))
  seqs <- ape::read.FASTA(test_path("test.fasta"))
  seq.info <- parse.headers(
    names(seqs), var.names=c('accn', 'coldate', 'subtype'),
    var.transformations=c(as.character, as.Date, as.factor))
  seq.info$colyear <- data.table::year(seq.info$coldate)
  which.new <- (seq.info$colyear >= 2012)
  obj <- read.edges(edge.info, seq.info, which.new)
  
  result <- component.cluster(obj, dist.thresh=0.06)
  expect_true(is.data.table(result))
  
  expected <- c(2, 1, 1, 2)  # (1,2), 3, 4, (5,6)
  expect_equal(result$Size, expected)
  expected <- c(1, 1, 0, 1)
  expect_equal(result$Growth, expected)
  expected <- list(
    c(as.Date("2008-11-04"), as.Date("2009-04-28")),
    as.Date("2011-02-01"),
    as.Date("2008-09-30"),
    c(as.Date("2010-03-02"), as.Date("2008-05-20"))
  )
  expect_equal(result$coldate, expected)
  
  # relaxing threshold to merge nodes 3 and 4
  result <- component.cluster(obj, dist.thresh=0.085)
  expect_equal(result$Size, c(2,2,2))
  expect_equal(result$Growth, c(1,1,1))
  result <- component.cluster(obj, dist.thresh=0.086)
  expect_equal(result$Size, c(4,2))
  expect_equal(result$Growth, c(2,1))
  
  # tightening threshold causes new case to drop out
  result <- component.cluster(obj, dist.thresh=0.05)
  expect_equal(result$Size, c(2, 1, 1, 2))
  expect_equal(result$Growth, c(1, 0, 0, 1))
  
  result <- component.cluster(obj, dist.thresh=0.04)
  expect_equal(result$Size, c(1, 1, 1, 1, 1, 1))
  expect_equal(result$Growth, c(0, 1, 0, 0, 1, 0))
})

