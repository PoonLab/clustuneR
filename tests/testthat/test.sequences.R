setwd("../../")

source("R/sequence.setup.R")
library(testthat)
library(ape)

load("data/seq.info_ex.RData")
load("data/seq-phylo_ex.RData")

test_that("Sequence Headers will be overwritten by full headers if named by User", {
  expect_warning(x <- pull.headers(seqs.full,var.names = c("Header", "CollectionDate", "Subtype"),
                              var.transformations =list(as.character, as.Date, as.factor)),
                 "'Header' is contained within var.names, this will be overwritten")
  expect_true(all(x$Header%in%names(seqs.full)))
})



