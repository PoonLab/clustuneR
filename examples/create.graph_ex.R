load("data/seq-phylo_ex.RData")
load("data/seq.info_ex.RData")

edge.info <- ape::dist.dna(seqs.full, pairwise.deletion = T, as.matrix = T, model = "TN93", )

new.year <- max(seq.info$CollectionDate) - 365
which.new <- which(seq.info$CollectionDate > new.year)
g <- create.graph(seq.info, edge.info, which.new)
