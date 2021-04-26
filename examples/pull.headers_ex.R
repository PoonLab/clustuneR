load("data/seq-phylo_ex.RData")

seq.info <- pull.headers(seqs.full,var.names = c("ID", "CollectionDate", "Subtype"),
                         var.transformations =list(as.character, as.Date, as.factor))
