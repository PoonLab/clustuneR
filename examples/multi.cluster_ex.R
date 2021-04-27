param.list.boot <- lapply(c(0.005,0.01), function(x){
  list("t"=extended.tree.ex, "branch.thresh"=x, "boot.thresh"=0.7)
})

step.cluster.range <- multi.cluster(step.cluster, param.list.boot, mc.cores = 1)
