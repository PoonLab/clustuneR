setwd("../../")

source("R/analysis.R")
library(testthat)
library(ape)

load("data/clusters_ex.RData")
load("data/fit.results_ex.RData")

predictive.models <- list(
  "NullModel" = function(x){
    glm(Growth~Size, data=x, family="poisson")
  },
  "TimeModel" = function(x){
    glm(Growth~Size+CollectionDate, data=x, family="poisson")
  }
)

predictor.transformations <- list(
  "CollectionDate" = function(x){mean(x)}
)


test_that("Cluster analysis can be performed on a single set of component.clusters", {
  expect_warning(expect_error(fit.analysis(component.cluster.set,
                              predictor.transformations = predictor.transformations,
                              predictive.models=predictive.models),NA),
                 "No range ID, by default this will be set to 0 for all sets")
})

test_that("Cluster analysis can be performed with no provided model", {
  cluster.data <- step.cluster.data
  expect_error(fit.analysis(cluster.data),NA)
})

