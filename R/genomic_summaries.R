#' Estimates minor allele frequency (MAF)
#'
#' @param M The additive \eqn{n \times p} matrix coded with 0, 1, 2, or \code{NA}
#' (default = \code{NULL}).
#'
#' @return A vector with the MAF of each molecular marker
#'
#' @keywords internal

maf <- function(M = NULL){

  # Get the frequency of the markers.
  maf <- colMeans(M, na.rm = TRUE)/2

  # Identify the MAF.
  maf <- apply(cbind(maf, 1 - maf), 1, min)

  return(maf)
}

#' Estimates observed and expected heterozygosity
#'
#' @param M The additive \eqn{n \times p} matrix coded with 0, 1, 2 coding (default = \code{NULL}).
#'
#' @return A list with vectors containing the observed (ho) and expected (he) heterozygosity.
#'
#' @keywords internal

heterozygosity <- function(M = NULL){

  # Get q.
  q <- maf(M)

  # Get p.
  p <- 1 - q

  # Get the expected heterozygosity.
  he <- 2 * p * q

  # Get the obseved heterozygosity.
  ho <- colMeans(M == 1, na.rm = TRUE)

  return(data.frame(ho = ho, he = he))
}

#' Estimates call rate
#'
#' @param M The additive \eqn{n \times p} matrix with any coding (default = \code{NULL}).
#' @param margin A character indicating the margin for call rate calculations.
#' Options are: \code{row} and \code{col} (default = \code{row}).
#'
#' @return A vector containing the call rate.
#'
#' @keywords internal

callrate <- function(M = NULL, margin = c("row", "col")){

  # Collect input.
  margin <- match.arg(margin)

  # CR by row.
  if (margin == "row") cr <- 100 - rowSums(is.na(M))/ncol(M) * 100

  # CR by col.
  if (margin == "col") cr <- 100 - colSums(is.na(M))/nrow(M) * 100

  return(cr)

}

#' Estimates the population level inbreeding (Fis) by marker
#'
#' @param M The additive \eqn{n \times p} matrix with any coding (default = \code{NULL}).
#' @param margin A character indicating the margin for call rate calculations.
#' Options are: \code{row} and \code{col} (default = \code{col}).
#'
#' @return A vector containing the Fis for the markers.
#'
#' @keywords internal

Fis <- function(M = NULL, margin = c("col", "row")){

  # Collect input.
  margin <- match.arg(margin)

  # H by row.
  if (margin == "col") H <- heterozygosity(M = M)

  # H by col.
  if (margin == "row") H <- heterozygosity(M = t(M))

  # Calculate Fis.
  Fis <- ifelse(test = H[, "he"] == 0,
                yes = 0,
                no = 1 - (H[,"ho"] / H[,"he"]))

  return(Fis)

}
