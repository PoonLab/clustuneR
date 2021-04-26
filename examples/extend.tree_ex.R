load("data/seq-phylo_ex.RData")
load("data/seq.info_ex.RData")

t <- extend.tree(tree.old, seq.info, full.align = seqs.full, log.file = "data/IQTREE_log_ex.txt")
