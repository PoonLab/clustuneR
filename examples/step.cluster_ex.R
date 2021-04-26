load("data/tree_ex.RData")

step.cluster.set <- step.cluster(t, branch.thresh = 0.007, boot.thresh = 0.90)
