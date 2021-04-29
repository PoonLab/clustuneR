#' Get seq.info from headers
#'
#' Translates a set of sequence headers into a data.frame object for data input.
#' For well formatted alignments, this may be an easy way to obtain meta-data.
#' 
#' @param seqs: An inputted alignment using ape's sequence handling
#' @param var.names: The names of the variables represented in each header.
#' @param var.transformations: A list of transformation functions (such as as.Date())
#' these transform each row into it's proper type. by default, each type is character.
#' @param sep: The separator character that splits up the headers into variables
#' @return A data.table object containing sequence meta data paired to headers. 
#' See seq.info.ex, for more details and an example
#' @export
#' @example examples/pull.headers_ex.R
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
