step.cluster.set <- step.cluster(extended.tree.ex, branch.thresh = 0.03, boot.thresh = 0)

step.cluster.set[which.max(Size), ]
