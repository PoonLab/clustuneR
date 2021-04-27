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

cluster.data <- cluster.ex
cluster.data[,"RangeID":=0]

res <- fit.analysis(cluster.data,
                    predictor.transformations = predictor.transformations,
                    predictive.models = predictive.models)

aics <- get.AIC(res)
aicdiff <- aics$TimeModelAIC-aics$NullModelAIC
which.min(aicdiff)
