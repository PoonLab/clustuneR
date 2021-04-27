cluster.data <- cluster.ex
cluster.data[,"RangeID":=0]

fit.result <- fit.analysis(cluster.data)

mod.performance <- fit.result$NullModel
