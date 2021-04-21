#' Multiple clusters from a parameter set
#'
#' Runs a given clustering method over a range of parameters values.
#' If growth is attached This result can be fed into a fit.analysis
#'
#' @param param.list: A named list of parameter sets. Each must correspond to the clustering method used.
#' @param rangeID: If several different parameter ranges are used, the rangeID can identify them.
#' @param mc.cores: A parallel option
#' @return A larger data.table with parameter sets noted
multi.cluster <- function(cluster.method, param.list, mc.cores = 1, rangeID = 0) {

  # Cluster method loop
  cluster.range <- parallel::mclapply(1:length(param.list), function(i) {
    x <- param.list[[i]]
    x$setID <- i
    do.call(cluster.method, x)
  }, mc.cores = mc.cores)

  cluster.range <- dplyr::bind_rows(cluster.range)
  suppressWarnings(cluster.range[, "RangeID" := rangeID])

  return(cluster.range)
}

#' Predictive analysis on a range of cluster sets
#'
#' Fits a predictive model of some outcome (by default, cluster growth) to sets of cluster data.
#' This fit measurement is recorded for each use of the predictive model on a given cluister set
#'
#' @param cluster.data: Inputted set(s) of cluster data. May or may not be sorted into ranges
#' @param mc.cores: A parallel option
#' @param predictor.transformations: A named list of transformation functions for each predictor variable.
#' This name should correspond to a column from the cluster.data, which will be taken as input for the function
#' @param predictive.models: A named list of functions, each of which applies a model to inputted data (x). See default null for example.
#' @return A data.table of analysis results. Several important summary values such as null and full AIC are proposed here
fit.analysis <- function(cluster.data, mc.cores = 1, predictor.transformations = list(),
                         predictive.models = list(
                           "NullModel" = function(x){
                             glm(Size~Growth, data=x, family="poisson")
                             })) {
  # Check inputs
  predictors <- names(predictor.transformations)
  mod.names <- names(predictive.models)
  setIDs <- unique(cluster.data[, SetID])
  if (!all((predictors) %in% colnames(cluster.data))) {
    stop("Predictors referenced in transform step are not in the range of cluster data")
  }
  if (!("Growth" %in% colnames(cluster.data))) {
    warning("No Growth information with clusters. This will default to 0 for all clusters.
            Some other outcome should be specified in your predictive models.")
    cluster.data[, "Growth" := 0]
  }
  if (!("RangeID" %in% colnames(cluster.data))) {
    warning("No range ID, by default this will be set to 0 for all sets")
    cluster.data[, "RangeID" := 0]
  }

  # Transform cluster data for modelling based on inputs
  model.data <- cluster.data[, c("Header", "Size", "Growth", "SetID", "RangeID")]

  if(!is.null(predictors)) {
    model.data[, (predictors) := lapply(predictors, function(x) {
      sapply(cluster.data[, get(x)], function(z) {
        (predictor.transformations[[x]])(z)
      })
    })]
  }


  # Obtain fit data for each cluster set
  cluster.analysis <- dplyr::bind_rows(
    parallel::mclapply(setIDs, function(id) {
      DT <- model.data[SetID == id, ]

      res <- data.table::data.table("SetID" = DT[1, SetID], "RangeID" = DT[1, RangeID])
      res[, (mod.names) := lapply(predictive.models, function(pmod){suppressWarnings(list(pmod(DT)))})]
      return(res)
    }, mc.cores = mc.cores)
  )

  return(cluster.analysis)
}

#'Get AIC differences from an analysis
#'
#'Takes a cluster.analysis and extracts AIC values from columns containing models
#'
#'@param cluster.analysis: A data.table from some predictive growth model analysis
#'@return The AIC data for all columns containing fit objects
get.AIC <- function(cluster.analysis){

  #Identify models
  which.models <- sapply(cluster.analysis[1,], function(x){any(attr(x[[1]], "class")%in%c("lm", "glm"))})
  which.models <- which(which.models)
  if(length(which.models)==0) {
    stop("No models in the data set provided")
  }
  model.fits <- cluster.analysis[,.SD, .SDcols = which.models]

  #Return specifically aic as a new data.table
  newnms <- sapply(names(which.models), function(nm){paste0(nm,"AIC")})
  aic.analysis <- data.table::data.table()
  aic.analysis[,(newnms) := lapply(names(which.models), function(nm){
    sapply(model.fits[,get(nm)], function(x){x$aic})
  })]

  return(aic.analysis)
}
