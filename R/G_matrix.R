#' Obtains the genomic matrix from SNP data for additive or dominant relationships
#'
#' Generates the genomic numerator relationship matrix for
#' additive (VanRaden or Yang) or dominant (Su or Vitezica) relationships.
#' Matrix provided \eqn{\boldsymbol{M}} is of the form \eqn{n \times p}, with \eqn{n} individuals and \eqn{p} markers.
#' Individual and
#' marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as 0, 1, 2 (integers or decimal numbers).
#' Missing values, if present, need to be specified.
#'
#' @param M A matrix with SNP data of form \eqn{n \times p}, with \eqn{n} individuals and \eqn{p} markers.
#' Individual and marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as 0, 1, 2 (integers or decimal numbers) (default = \code{NULL}).
#' @param method The method considered for calculation of genomic matrix.
#' Options are: \code{"VanRaden"} and \code{"Yang"} for additive matrix,
#' and \code{"Su"} and \code{"Vitezica"} for dominant matrix (default =  \code{"VanRaden"}).
#' @param na.string A character that is interpreted as missing values (default = \code{"NA"}).
#' @param sparseform If \code{TRUE} it generates a matrix in sparse form to be used directly in \pkg{asreml} with
#' required attributes (default = \code{FALSE}).
#' @param digits Set up the number of digits used to round the output matrix (default = 8).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @details
#' Note: If data is provided with missing values, it will process calculations
#' of relationships on pairwise non-missing data.
#'
#' It uses function \code{Gmatrix} for calculations
#' from R package  \pkg{AGHmatrix} (Amadeu \emph{et al.} 2019).
#'
#' @return A list with one of these two elements:
#' \itemize{
#' \item \code{G}: the \eqn{\boldsymbol{G}} matrix in full form (only if \code{sparseform = FALSE}).
#' \item \code{G.sparse}: the \eqn{\boldsymbol{G}} matrix in sparse form (only if \code{sparseform = TRUE}).
#' }
#'
#' @references
#' Amadeu, R.R., Cellon, C., Olmstead, J.W., Garcia, A.A.F, Resende, M.F.R. and P.R. Munoz. 2016.
#' AGHmatrix: R package to construct relationship matrices for autotetraploid and diploid species:
#' A blueberry example. The Plant Genome 9(3). doi: 10.3835/plantgenome2016.01.0009
#'
#' @export
#'
#' @examples
#' # Example: Requesting a full matrix by VanRanden.
#'
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#' G[1:5, 1:5]
#'
#' \donttest{
#' # Example: Requesting a sparse form by VanRanden.
#'
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden", sparseform = TRUE)$G.sparse
#' head(G)
#' head(attr(G, "rowNames"))
#' }
#'

G.matrix <- function(M = NULL, method = "VanRaden", na.string = "NA",
                     sparseform = FALSE, digits = 8, message = TRUE){

  use.dll <- FALSE
  # Check if the class of M is matrix
  if (is.null(M) || !inherits(M, "matrix")) {
    stop("M should be a valid object of class matrix.")
  }
  M <- as.matrix(M)
  n <- dim(M)[1]
  p <- dim(M)[2]
  if (p < n) {warning('There are more individuals than markers in this data!')}

  # Check the rownames/colnames of M
  if(is.null(colnames(M))) {
    stop('Marker names not assigned to columns of matrix M.')
  }
  if(is.null(rownames(M))) {
    stop('Individual names not assigned to rows of matrix M.')
  }

  # Call libraries to obtain G matrix
  # Using AGHmatrix
  if (!use.dll) {
    #options(echo=FALSE,message=FALSE)
    if (message){
      GHAT <- Gmatrix(M, method=method, missingValue=na.string, maf=0)
    } else {
      GHAT <- silent_(Gmatrix(M, method=method, missingValue=na.string, maf=0))
    }

    #options(echo=TRUE)
  # Using our own DLL
  } else {
    # # Preparing things for the DLL
    # for (i in 1:p) { M[,i] <- as.double(M[,i]) }  # Coersing M to be double precision
    # M <- as.matrix(M)
    # if (method == 'VanRaden') {mm = 2}
    # if (method == 'Yang') {mm = 1}
    # if (method == 'Su') {stop('Method not available in DLL.')}
    # if (method == 'Vitizica') {stop('Method not available in DLL.')}
    # # dyn.load("ghat4.dll", type='Fortran')
    # dyn.load("src/ghat4.dll", type='Fortran')  # Not sure this is the best way!
    # GHAT <- base::.Fortran("ghat", molM=M, ghatM=matrix(as.double(0),n,n),
    #                        mm=as.integer(2), n=as.integer(n), m=as.integer(p))$ghatM
    # #dyn.unload("ghat4.dll")
  }

  # round GHAT
  GHAT <- round(GHAT, digits)

  # Making sure GHAT has the rownames and colnames
  #GHAT <- as.matrix(GHAT)
  if (length(rownames(M)) == 0) {
    rownames(GHAT) <- as.character(1:n)
    colnames(GHAT) <- as.character(1:n)
  } else {
    rownames(GHAT) <- rownames(M)
    colnames(GHAT) <- rownames(M)
  }

  if (isTRUE(sparseform)) {
    GHAT.sparse <- full2sparse(GHAT)
    return(list(G.sparse = GHAT.sparse))
    } else{
    return(list(G = GHAT))
    }


}
