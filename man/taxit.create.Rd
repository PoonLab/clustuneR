% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree-setup.R
\name{taxit.create}
\alias{taxit.create}
\title{Create a refpkg}
\usage{
taxit.create(t, seqs.full, stats.json, locus = "LOCUS")
}
\arguments{
\item{t:}{The tree (made on a subset of the full alignment)}

\item{seqs.full:}{The full alignment. Including sequences excluded from the tree.}

\item{stats.json:}{Path to the full alignment file (new + old seqs)}

\item{locus:}{Extra information required for the summary json}
}
\value{
Path to a temporary refpkg directory
}
\description{
A wrapper for the taxit.create function used by pplacer. This will generate a
temporary refpkg file based on passed tree and alignment data summary json.
}
