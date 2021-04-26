#' Get seq.info from headers
#'
#' Translates a set of sequence headers into a data.frame object for data input.
#' NOTE: This must contain, at minimum, a set of unique sequence id's (labelled ID) and one other variable
#'
#' @param seqs: An inputted alignment using ape's sequence handling
#' @param var.names: The names of the variables represented in each header. This must contain "ID".
#' @param var.transformations: A list of transformation functions (such as as.character())
#' these transform each row into it's proper type. by default, each type is set to character.
#' @param sep: The separator character that splits upthe headers in the fasta file
#' @return A data.table object containing the information associated with each sequence.
#' @export
#' @examples
#' load("data/seq-phylo_ex.RData")
#' seq.info <- pull.headers(seqs.full,var.names = c("ID", "CollectionDate", "Subtype"),var.transformations =list(as.character, as.Date, as.factor))
pull.headers <- function(seqs, var.names, var.transformations = list(), sep = "_") {

  # Checking Inputs
  if (length(var.names) != length(unique(var.names))) {
    stop("var.names may not contain repeats")
  }
  if (length(var.transformations) == 0) {
    var.transformations <- lapply(1:length(var.names), function(x) {
      as.character
    })
  } else {
    if (length(var.names) != length(var.transformations)) {
      stop("var.names and var.transformations must be equal lengths")
      return(NULL)
    }
  }
  if ("Header" %in% var.names) {
    warning("'Header' is contained within var.names, this will be overwritten")
  }


  # Split and transform data from headers
  split.headers <- sapply(names(seqs), function(x) {
    strsplit(x, sep)[[1]]
  })


  seq.info <- lapply(1:nrow(split.headers), function(i) {
    x <- unname(split.headers[i, ])
    x <- var.transformations[[i]](x)
    data.table::data.table(x)
  })

  seq.info <- dplyr::bind_cols(seq.info)
  colnames(seq.info) <- var.names

  seq.info[, "Header" := names(seqs)]

  return(seq.info)
}

#' Obtain "old" sequences
#'
#' Get a subset of sequences from an ape seq object based on seq.info data which is not labeled new
#'
#' @param seqs: A full alignment
#' @param seq.info: A set of seq.info. If provided, may be used to ascertain old sequences from new
#' A tree is generally only to be built from old sequences
#' @return Filtered sequences, an ape seq object
get.old.seqs <- function(seqs, seq.info) {

  # Check new sequences, filter if given
  if ("New" %in% colnames(seq.info)) {
    stop("No new info given in seq.info")
  } else {
    seqs <- seqs[which(names(seqs) %in% seq.info[!(New), Header])]
  }

  return(seqs)
}
