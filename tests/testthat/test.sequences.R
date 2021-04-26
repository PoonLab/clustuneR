setwd("../../")

source("R/sequence.setup.R")
library(testthat)
library(ape)
library(phangorn)

load("data/seq.info_ex.RData")
load("data/seq-phylo_ex.RData")

test_that("Sequence Headers can be pulled from sequence data", {
  expect_error(pull.headers(seqs.full,var.names = c("ID", "CollectionDate", "Subtype"),
                             var.transformations =list(as.character, as.Date, as.factor)), NA)
})



