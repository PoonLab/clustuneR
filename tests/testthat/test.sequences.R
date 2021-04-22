source("R/sequence.setup.R")
library(testthat)

load("data/seq.info_ex.RData")
load("data/seq-phylo_ex.RData")

test_that("Sequence Headers can be pulled from sequence data", {
  expect_error(pull.headers(seqs.full,var.names = c("ID", "CollectionDate", "Subtype"),
                             var.transformations =list(as.character, as.Date, as.factor)), NA)
})

test_that("Sequences can be identified as new", {
  new.year <- max(seq.info$CollectionDate) - 365
  which.new <- which(seq.info$CollectionDate > new.year)
  expect_error(annotate.new(seq.info,which.new), NA)
})



