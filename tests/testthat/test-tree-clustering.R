# prepare test fixture
seqs <- ape::read.FASTA(test_path("test.fasta"))
seq.info <- parse.headers(
  names(seqs), var.names=c("accn", "coldate", "subtype"),
  var.transformations = c(as.character, as.Date, as.factor))
phy <- ape::read.tree(test_path("test-oldseq.fasta.treefile"))
phy <- import.tree(phy, seq.info = seq.info)
log.file <- test_path("test-oldseq.fasta.log")
obj <- extend.tree(phy, seqs, log.file)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
