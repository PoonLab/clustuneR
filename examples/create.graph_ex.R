edge.info <- ape::dist.dna(alignment.ex, pairwise.deletion = T, as.matrix = T, model = "TN93", )

new.year <- max(seq.info.ex$CollectionDate) - 365
which.new <- which(seq.info.ex$CollectionDate > new.year)
g <- create.graph(seq.info.ex, edge.info, which.new)

component.cluster(g)
