#' Check any kinship matrix K against phenotypic data
#'
#' Assesses a given kinship matrix against the provided phenotypic data to determine
#' if all genotypes are in the kinship matrix or not. It also reports which individuals
#' match or are missing from one set or another.
#' If requested, a reduced kinship matrix is generated that has only the matched individuals.
#' The input kinship matrix can be a pedigree-based relationship matrix \eqn{\boldsymbol{A}},
#' a genomic-based relationship matrix \eqn{\boldsymbol{G}}, or a hybrid
#' relationship matrix \eqn{\boldsymbol{H}}.
#' Individual names should be assigned to \code{rownames} and \code{colnames} of input matrix.
#'
#' @param K Input of a kinship matrix in full form (\eqn{n \times n}) (default = \code{NULL});
#' @param pheno.data A data fame with the phenotypic data to assess (for \code{n} individuals)
#' (default = \code{NULL}).
#' @param indiv The string for the column name for genotypes/individuals in the phenotypic data (default = \code{NULL}).
#' @param clean If \code{TRUE}, generates a new clean kinship matrix containing only the matched
#' phenotyped individuals (default = \code{FALSE}).
#' @param ord If \code{TRUE}, it will order the kinship matrix as in the phenotypic data, which is
#' recommended for some downstream genomic analyses (default = \code{TRUE}).
#' @param mism If \code{TRUE}, it generates data frames with matched and mismatched individual's names
#' from the kinship matrix and the phenotypic data (default = \code{FALSE}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with the following elements:
#'  \itemize{
#' \item{\code{mismatchesK}: a vector containing the names of the individuals from the provided kinship matrix
#' that \emph{mismatch} with the phenotypic data.}
#' \item{\code{matchesK}: a vector containing the names of the individuals from the provided kinship matrix
#' that \emph{match} with the phenotypic data.}
#' \item{\code{mismatchesP}: a vector containing the names of phenotyped individuals
#' that \emph{mismatch} with those from the kinship matrix.}
#' \item{\code{matchesP}: a vector containing the names of phenotyped individuals
#' that \emph{match} with those from the kinship matrix.}
#' \item{\code{Kclean}: a clean kinship matrix containing only the matched phenotyped individuals.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get G matrix.
#' G <- G.matrix(M = geno.pine655, method = "VanRaden", na.string = "-9", sparseform = FALSE)$G
#' dim(G)
#'
#' # Match G and the phenotype.
#' check <-
#'  match.kinship2pheno(
#'   K = G, pheno.data = pheno.pine,
#'   indiv = "Genotype",
#'   clean = TRUE, mism = TRUE)
#' ls(check)
#' length(check$matchesK)
#' length(check$mismatchesK)
#' length(check$matchesP)
#' length(check$mismatchesP)
#' dim(check$Kclean)
#'}


match.kinship2pheno <- function(K = NULL,
                                pheno.data = NULL, indiv = NULL,
                                clean = FALSE, ord = TRUE, mism = FALSE,
                                message = TRUE){

  # Check if the class of K is matrix
  if (is.null(K) || !inherits(K, "matrix")) {
    stop("K should be a valid object of class matrix")
  }
  if (is.null(rownames(K))){
    stop("Individual names not assigned to rows of matrix K")
  }
  if (is.null(colnames(K))){
    stop("Individual names not assigned to columns of matrix K")
  }
  if ((identical(rownames(K), colnames(K))) == FALSE){
    stop("Rownames and colnames of matrix K do not match.")
  }

  # Checks
  kinship_names <- rownames(K)
  pheno_names <- pheno.data[[indiv]]

  if( all(kinship_names %in% pheno_names) == TRUE){
    if (message){
      message("All individuals within the kinship matrix match the phenotyped individuals.")
    }
    mismatchesK <- NULL
    matchesK <- which(kinship_names %in% pheno_names == TRUE)
  } else {
    mismatchesK <- which(kinship_names %in% pheno_names != TRUE)
    if (message){
      message("Kinship matrix contains ", length(mismatchesK),
              " individuals that DO NOT match the phenotyped individuals.")
    }
    matchesK <- which(kinship_names %in% pheno_names == TRUE)
    if (message){
      message("Kinship matrix contains ", length(matchesK),
              " individuals that match the phenotyped individuals.")
    }
  }

  if( all(pheno_names %in% kinship_names) == TRUE){
    if (message){
      message("All phenotyped individuals match the individuals within the kinship matrix.")
    }
    mismatchesP <- NULL
    matchesP <- which(pheno_names %in% kinship_names == TRUE)
  } else {
    mismatchesP <- which(pheno_names %in% kinship_names != TRUE)
    if (message){
      message("Phenotypic data contains ", length(mismatchesP),
              " individuals that DO NOT match the kinship matrix individuals.")
    }
    matchesP <- which(pheno_names %in% kinship_names == TRUE)
    if (message){
    message("Phenotypic data contains ", length(matchesP),
            " individuals that match the kinship matrix individuals.")
    }
  }

  if(length(mismatchesK) != 0 & clean){
    if (message){
      message("Individuals within the kinship matrix that do not match those in the phenotypic data will be removed from this matrix.")
    }
    #Kclean <- K[matchesP, matchesP]
    Kclean <- K[matchesK, matchesK]
  } else {
    Kclean <- NULL
  }

  if(ord) {
    #list.order <- match(rownames(Kclean), pheno_names)
    list.order <- match(pheno_names, rownames(Kclean))
    Kclean <- Kclean[list.order,list.order]
  }
  if(!mism) {
    mismatchesK <- NULL; matchesK <- NULL
    mismatchesP <- NULL; matchesP <- NULL
  }

  return(list(mismatchesK=mismatchesK, matchesK=matchesK,
              mismatchesP=mismatchesP, matchesP=matchesP,
              Kclean=Kclean))

}
