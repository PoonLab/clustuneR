test_that("read.edges works", {
  edge.info <- read.csv(test_path("test.tn93.csv"))
  
  seqs <- ape::read.FASTA(test_path("test.fasta"))
  seq.info <- parse.headers(
    names(seqs), var.names=c('accn', 'coldate', 'subtype'),
    var.transformations=c(as.character, as.Date, as.factor))
  seq.info$colyear <- data.table::year(seq.info$coldate)
  
  idx <- match(unique(edge.info$ID1), seq.info$Header)
  edge.info <- edge.info[]
  
  which.new <- (seq.info$colyear >= 2012)
  expect_equal(sum(which.new), 3)
  
  obj <- read.edges(edge.info, seq.info, which.new, 
                    growth.resolution=minimum.retrospective.edge)
  
  # note, edge resolution reduces number of edges
  result <- obj$edge.info$ID1
  expect_equal(length(result), (6*5)/2 + 3)
  
  expected <- c(
    rep(1, 5),  # KU190000
    rep(2, 4),  # KU190006 
    rep(3, 3),  # KU190950 
    rep(4, 2),  # KU190868
    rep(5, 1),  # KU190389 
    5, # KU190389 + KU190613
    2, # KU190006 + KU190027
    3  # KU190950 + KU191003
    )
  expect_equal(result, expected)
  
  expected <- c(2:6, 3:6, 4:6, 5:6, 6, 7:9)
  expect_equal(obj$edge.info$ID2, expected)
})

