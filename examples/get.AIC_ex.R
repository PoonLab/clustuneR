load("data/fit.results_ex.RData")
DT <- get.AIC(fit.result)
aicdiff <- DT$TimeModelAIC-DT$NullModelAIC
which.min(aicdiff)