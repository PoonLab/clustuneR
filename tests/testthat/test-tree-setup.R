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
