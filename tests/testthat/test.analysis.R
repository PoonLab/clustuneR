setwd("../../")

source("R/analysis.R")
library(testthat)
library(ape)
library(phangorn)

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

test_that("Cluster analysis can be performed on a set of extreme step.clusters", {
  expect_error(fit.analysis(step.cluster.data,
                            predictor.transformations = predictor.transformations,
                            predictive.models=predictive.models),NA)
})

test_that("Cluster analysis can be performed on a set of extreme mono.pat.clusters which have no growth or variables", {
  expect_error(fit.analysis(mono.pat.cluster.data),NA)
})

test_that("Cluster analysis can be performed on a set of extreme component.clusters", {
  expect_error(fit.analysis(component.cluster.data,
                            predictor.transformations = predictor.transformations,
                            predictive.models=predictive.models),NA)
})

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

test_that("AIC can be pulled from analysis sets", {
  expect_error(get.AIC(res.step.boot), NA)
  expect_error(get.AIC(res.step.noboot), NA)
})

