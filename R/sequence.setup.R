#' Convenience function to parse metadata from sequence headers
#'
#' Translates a set of sequence headers into a data.frame object for data input.
#' Requires that every sequence header contains the same number and type of 
#' values, separated by a reserved delimiter character, e.g., underscore "_"
#'  
#' @param headers: character, vector of sequence names such as obtained by 
#' names(x) where x is an ape::DNAbin object.
#' @param var.names: character, The names of the variables represented in each header.
#' @param var.transformations: list, A list of transformation functions (such as as.Date())
#' these transform each row into it's proper type.  By default, each value is taken 
#' literally (i.e., as.character).
#' @param sep: character, The separator character that splits up the headers into 
#' values.
#' @return A data.table object containing sequence meta data paired to headers. 
#' See seq.info.ex, for more details and an example
#' @export
#' @example examples/parse.headers_ex.R
parse.headers <- function(headers, var.names=c(), var.transformations = list(), sep = "_") {

  # Split and transform data from headers
  split.headers <- strsplit(headers, sep)
  raw.mx <- do.call(rbind, split.headers)
  
  # Checking Inputs
  if (length(var.names) == 0) {
    warning("Warning, you did not specify any variable names. Using default V1, V2, etc.")
    var.names <- paste("V", 1:ncol(raw.mx), sep="")
  }
  if (length(var.names) != length(unique(var.names))) {
    stop("var.names may not contain repeats")
  }
  if (length(var.transformations) == 0) {
    var.transformations <- sapply(1:length(var.names), function(x) {
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
  
  seq.info <- data.frame(Header=headers)  # store original values
  for (i in 1:ncol(raw.mx)) {
    seq.info[var.names[i]] <- var.transformations[[i]](raw.mx[,i])
  }
  
  return(seq.info)
}
