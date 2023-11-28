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
})

test_that("annotate.paths works", {
  phy <- ape::read.tree(text=nwk)
  phy <- phangorn::midpoint(phy)
  phy <- ape::multi2di((phy))
  phy$node.info <- annotate.nodes(phy)
  
  result <- annotate.paths(phy)
  expect_type(result, 'list')
})
