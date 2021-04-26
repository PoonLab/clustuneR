source("R/tree.clustering.R")
load("data/tree_ex.RData")

param.list.boot <- lapply(c(0.005,0.01), function(x){
  list("t"=t, "branch.thresh"=x, "boot.thresh"=0.9)
})

step.cluster.data.boot <- multi.cluster(step.cluster, param.list.boot, mc.cores = 1)
