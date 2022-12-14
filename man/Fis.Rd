% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomic_summaries.R
\name{Fis}
\alias{Fis}
\title{Estimates the population level inbreeding (Fis) by marker}
\usage{
Fis(M = NULL, margin = c("col", "row"))
}
\arguments{
\item{M}{The additive \eqn{n \times p} matrix with any coding (default = \code{NULL}).}

\item{margin}{A character indicating the margin for call rate calculations.
Options are: \code{row} and \code{col} (default = \code{col}).}
}
\value{
A vector containing the Fis for the markers.
}
\description{
Estimates the population level inbreeding (Fis) by marker
}
\keyword{internal}
