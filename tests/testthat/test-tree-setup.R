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
  
  result <- annotate.nodes(phy)
  expect_true(is.data.table(result))
  
  expected <- c(0.9, 0.99, 1.0, 1.0, 1.0)
  expect_equal(sort(result$Bootstrap), expected)  
})

test_that("annotate.paths works", {

})
