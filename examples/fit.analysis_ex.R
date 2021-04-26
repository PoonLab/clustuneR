load("data/clusters_ex.RData")

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

fit.result <- fit.analysis(component.cluster.data,
                           predictor.transformations = predictor.transformations,
                           predictive.models=predictive.models)
