nwk <- "(KU190000.1_2008-11-04_subA:0.0207974896,KU190006.1_2009-04-28_subA:0.0376395280,\
((KU190950.1_2011-02-01_subC:0.0392260558,KU190868.1_2008-09-30_subC:0.0749029434)90:\
0.0257377758,(KU190389.1_2010-03-02_subB:0.0327358095,KU190083.1_2008-05-20_subB:\
0.0233397233)100:0.0768084940)99:0.0360960580);"

test_that("import.tree works", {
  phy <- ape::read.tree(text=nwk)
  
  # call without passing seq.info
  expect_warning(result <- import.tree(phy))
  
  expect_true(is.rooted(result))
  expect_true(is.binary(result))
  expect_equal(6, Ntip(result))
  expect_equal(5, Nnode(result))
  })

test_that("annotate.nodes works", {
  phy <- ape::read.tree(text=nwk)
  phy <- phangorn::midpoint(phy)
  phy <- ape::multi2di((phy))
  
  # nodes are numbered 1..n for n tips, and n+1..n+m for m internal nodes
  # root is n+1, and rest of numbering by preorder traversal
  dt <- annotate.nodes(phy)
  expect_true(is.data.table(dt))
  
  # check that bootstrap values are rescaled to (0,1)
  expected <- c(0.9, 0.99, 1.0, 1.0, 1.0, rep(1.0, Ntip(phy)))
  expect_equal(sort(dt$Bootstrap), expected)
  
  # check numbers of descendants per node (does not count self internal)
  expected <- sort(c(10, 6, 2, 2, 2, rep(1, 6)))
  result <- sort(sapply(dt$Descendants, length))
  expect_equal(result, expected)
  
  # manually calculated patristic distances
  expect_equal(max(dt$max.patristic.dist), 0.210185, tolerance=1e-6)
  expected <- min(dt$max.patristic.dist[dt$max.patristic.dist>0])
  expect_equal(expected, 0.05607553, tolerance=1e-6)
  
  # cherries have only one patristic distance
  result <- sort(c(0.074902943 + 0.039226056, 0.037639528 + 0.020797490, 
              0.032735809 + 0.023339723))
  expected <- sort(dt$mean.patristic.dist)[1:3]
  expect_equal(result, expected, tolerance=1e-6)
  
  expect_equal(dt$NodeID, 1:11)
  
  expected <- c(0.020797490, 0.037639528, 0.039226056, 0.074902943,
                0.032735809, 0.023339723, NA, 0.036096058, 0.004451792,
                0.025737776, 0.072356702)
  expect_equal(dt$BranchLength, expected)
  
  # check descendants
  expected <- list(1, 2, 3, 4, 5, 6, c(1:6, 8:11), c(1,2), c(1:4, 8, 10), 
                   c(3,4), c(5,6))
  for (i in 1:nrow(dt)) {
    expect_setequal(dt$Descendants[[i]], expected[[i]])  
  }
  
  # check paths to root
  expected <- list(
    c(1, 8, 9, 7), c(2, 8, 9, 7), c(3, 10, 9, 7), c(4, 10, 9, 7), c(5, 11, 7), 
    c(6, 11, 7), c(7), c(8, 9, 7), c(9, 7), c(10, 9, 7), c(11, 7)
    )
  expect_equal(dt$Paths, expected)
})


test_that("translate.log works", {
  json <- translate.log(test_path("test-oldseq.fasta.log"))
  expect_true(class(json) == 'json')
  result <- jsonlite::fromJSON(json)
  expect_true(class(result) == 'list')
  
  # TODO: some of these values are constant, not worth checking
  expect_true(result$empirical_frequencies)
  expect_equal(result$datatype, "DNA")
  expect_equal(result$subs_model, "GTR")
  expect_equal(result$program, "IQ-TREE")
  
  expected <- list("ac"=2.43167, "gt"=1.00000, "at"=1.00000, 
                   "ag"=12.56081, "cg"=2.43167, "ct"=12.56081)
  expect_equal(result$subs_rates, expected, tolerance=0.0001)
})

test_that("annotate.growth works", {
  seqs <- ape::read.FASTA(test_path("test.fasta"))
  seq.info <- parse.headers(
    names(seqs), var.names=c("accn", "coldate", "subtype"),
    var.transformations = c(as.character, as.Date, as.factor))
  
  phy <- ape::read.tree(test_path("test-oldseq.fasta.treefile"))
  phy <- import.tree(phy, seq.info = seq.info)
  
  # bypasses first part of extend.tree()
  ptrees <- ape::read.tree(test_path("growth.tre"))
  ptree <- ptrees[[1]]
  
  result <- annotate.growth(phy, ptree)
  expect_true(is.data.table(result))
  expect_equal(ncol(result), 6)
  
  # check contents of result against manually curated values
  expected <- rep("KU191003.1_2012-11-19_subC", 2)
  expect_equal(result$Header, expected)
  
  # KU190950.1_2011-02-01_subC, KU190868.1_2008-09-30_subC
  expected <- c(3, 4)
  expect_setequal(result$NeighbourDes, expected)
  
  expected <- c(0.0206413, 0.979359)
  expect_equal(sort(result$Bootstrap), expected, tolerance=1e-6)
  
  expected <- c(0.027342, 0.031319)
  expect_equal(sort(result$TermDistance), expected, tolerance=1e-5)
  
  expected <- c(0.031486, 0.0652853)
  expect_equal(sort(result$PendantDistance), expected, tolerance=1e-5)
  
  expected <- c(TRUE, TRUE)
  expect_equal(result$Terminal, expected)
  
  # check another placement tree
  ptree <- ptrees[[2]]
  result <- annotate.growth(phy, ptree)
  expect_equal(result$Header, "KU190027.1_2012-11-22_subA")
  expect_equal(result$NeighbourDes, list(2))
  expect_equal(result$Bootstrap, 1)
  expect_equal(result$PendantDistance, 0.0237015, tolerance=1e-6)
  expect_equal(result$Terminal, TRUE)
})

test_that("extend.tree works", {
  #seqs <- ape::read.FASTA(system.file("exdata/test.fasta", package="clustuneR"))
  seqs <- ape::read.FASTA(test_path("test.fasta"))
  expect_equal(length(seqs), 9)
  
  seq.info <- parse.headers(
    names(seqs), var.names=c("accn", "coldate", "subtype"),
    var.transformations = c(as.character, as.Date, as.factor))
  
  phy <- ape::read.tree(test_path("test-oldseq.fasta.treefile"))
  phy <- import.tree(phy, seq.info = seq.info)
  expect_true("node.info" %in% names(phy))
  
  log.file <- test_path("test-oldseq.fasta.log")
  result <- extend.tree(phy, seqs, log.file)
  expect_true("growth.info" %in% names(result))
  expect_true(is.data.table(result$growth.info))
  expect_equal(nrow(result$growth.info), 4)
  
  # test if tips below induced nodes (NeighbourDes) are collapsed to their MRCA
  expect_true(!is.element("NeighbourDes", names(result$growth.info)))
  
  # we need another example to really check if neighbour collapse is working
  seqs <- ape::read.FASTA(test_path("test2.fasta"))
  seq.info <- parse.headers(
    names(seqs), var.names=c("accn", "coldate", "subtype"),
    var.transformations = c(as.character, as.Date, as.factor))
  
  phy <- ape::read.tree(test_path("test2-old.fasta.treefile"))
  phy <- import.tree(phy, seq.info = seq.info)
  expect_equal(sum(phy$seq.info$New), 4)
  
  log.file <- test_path("test2-old.fasta.log")
  extended <- extend.tree(phy, seqs, log.file)
  expect_true("growth.info" %in% names(extended))
  
  ## needed to retrieve NeighbourDes values before calling collapse.neighbours()
  #stats.json <- translate.log(log.file)
  #refpkg <- taxit.create(phy, seqs, stats.json) 
  #ptrees <- run.pplacer_guppy(refpkg)
  ## copy growth.tre from tmp dir to test folder
  
  # use a test fixture
  ptrees <- read.tree(test_path("growth2.tre"))
  growth.info <- lapply(ptrees, function(x) annotate.growth(phy, x))
  expect_true("NeighbourNode" %in% names(extended$growth.info))
  result <- extended$growth.info$NeighbourNode
  expect_true(all(sapply(result, length)==1))
  expected <- c(6, 18, 17, 15, 14, 4, 15, 15)
  expect_equal(result, expected)
})
