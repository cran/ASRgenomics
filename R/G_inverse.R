#' Obtains the inverse of the genomic relationship matrix G
#'
#' @description
#' Generates the inverse of a genomic relationship matrix \eqn{\boldsymbol{G}} that is provided.
#' This input matrix should be of the full form (\eqn{n \times n})
#' with individual names assigned to
#' \code{rownames} and \code{colnames}. Several checks for the stability of the matrix are
#' presented based on the reciprocal conditional number.
#'
#' In case of an ill-conditioned matrix,
#' options of blending, bending or aligning before inverting are available.
#' These options will be deprecated (discontinued)
#' in future versions of \code{ASRgenomics}
#' as they can be better implemented in the function \code{G.tuneup()}.
#'
#' Based on procedures published by Nazarian and Gezan \emph{et al.} (2016).
#'
#' @param G Input of the symmetric genomic relationship matrix \eqn{\boldsymbol{G}} in
#' full form (\eqn{n \times n}), to obtain its inverse (default = \code{NULL}).
#' @param A Input of the pedigree relationship matrix \eqn{\boldsymbol{A}}
#' to perform blending or aligning in full form.
#' It should be of the same dimension as the \eqn{\boldsymbol{G}} matrix
#'  (\eqn{n \times n}) (default = \code{NULL}) (to be deprecated).
#' @param rcn.thr A threshold for identifying  the \eqn{\boldsymbol{G}} matrix as an ill-conditioned matrix.
#' Based on the reciprocal conditional number (default = \code{1e-12}).
#' @param blend If \code{TRUE} a "blending" with identity matrix \eqn{\boldsymbol{I}} or pedigree relationship matrix
#' \eqn{\boldsymbol{A}} (if provided) is performed (default = \code{FALSE}) (to be deprecated).
#' @param pblend If blending is requested this is the proportion of the identity matrix \eqn{\boldsymbol{I}} or
#' pedigree relationship matrix \eqn{\boldsymbol{A}} to blend for (default = \code{0.02}) (to be deprecated).
#' @param bend If \code{TRUE} a "bending" is performed by making the matrix near positive definite (default = \code{FALSE}) (to be deprecated).
#' @param eig.tol Defines relative positiveness (\emph{i.e.}, non-zero) of eigenvalues compared to the largest one.
#' It determines which threshold of eigenvalues will be treated as zero (default = \code{NULL}) (to be deprecated).
#' @param align If \code{TRUE} the genomic relationship matrix \eqn{\boldsymbol{G}} is aligned to the
#' pedigree relationship matrix \eqn{\boldsymbol{A}} (default = \code{FALSE}) (to be deprecated).
#' @param digits Set up the number of digits in used to round the output matrix (default = \code{8}).
#' @param sparseform If \code{TRUE} it generates an inverse matrix in sparse form to be used directly in \pkg{asreml} with
#' required attributes (default = \code{FALSE}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with three of the following elements:
#' \itemize{
#' \item{\code{Ginv}: the inverse of \eqn{\boldsymbol{G}} matrix in full form (only if \code{sparseform = FALSE})}.
#' \item{\code{Ginv.sparse}: the inverse of \eqn{\boldsymbol{G}} matrix in sparse form (only if \code{sparseform = TRUE})}.
#' \item{\code{status}: the status (\code{ill-conditioned} or \code{well-conditioned}) of
#' the inverse of \eqn{\boldsymbol{G}} matrix.}
#' \item{\code{rcn}: the reciprocal conditional number of the inverse of \eqn{\boldsymbol{G}} matrix.}
#' }
#'
#' @references
#' Nazarian A., Gezan S.A. 2016. GenoMatrix: A software package for pedigree-based
#' and genomic prediction analyses on complex traits. Journal of Heredity 107:372-379.
#'
#' @export
#'
#' @examples
#' # Example: An ill-conditioned matrix.
#'
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#' G[1:5, 1:5]
#'
#' # Get the inverse of G.
#' GINV <- G.inverse(G = G, bend = FALSE, blend = FALSE, align = FALSE)
#' GINV$Ginv[1:5, 1:5]
#' GINV$status
#'

G.inverse <- function(G = NULL, A = NULL, rcn.thr = 1e-12,
                      blend = FALSE, pblend = 0.02,
                      bend = FALSE, eig.tol = NULL, align = FALSE,
                      digits = 8, sparseform = FALSE, message = TRUE){

  # TODO check why eig.tol is not used in the function.

  # Deprecation traps ---------------------------------------------------------------------------

  if (!is.null(A) | blend | bend | align | !is.null(eig.tol)){
    warning("The arguments \'A', \'blend', \'pblend', \'bend', \'align', \'eig.tol' are still active",
            " but will be discontinued in future versions. Please use \'G.tuneup'",
            " to perform the respective actions.")
  }

  # Traps ---------------------------------------------------------------------------------------

  # Check if the class of G is matrix
  if (is.null(G) || !inherits(G, "matrix")) {
    stop("G should be a valid object of class matrix.")
  }
  # Check the rownames/colnames
  if (is.null(rownames(G))){
    stop("Individual names not assigned to rows of matrix G.")
  }
  if (is.null(colnames(G))){
    stop('Individual names not assigned to columns of matrix G.')
  }
  if ((identical(rownames(G), colnames(G))) == FALSE){
    stop("Rownames and colnames of matrix G do not match.")
  }

  # Checks for other input.
  if (pblend < 0 | pblend > 1) {
    stop("Specification of pblend must be between 0 and 1.")
  }

  if (rcn.thr <= 0) {
    stop("Value for rcn.thr must be positive.")
  }

  if(!is.null(eig.tol)){
    if (eig.tol <= 0) {
      stop("Value for eig.tol must be positive.")
    }
  }

  # Get dimensions.
  n <- dim(G)[1]
  p <- dim(G)[2]

  # Reciprocal Condition Number RCN
  rcn <- rcond(G)
  if (message){
    message('Reciprocal conditional number for original matrix is: ', rcn)
  }

  # Check the RCN default value and inverting the matrix
  if (rcn > rcn.thr){

    # Try to get the inverse of G.
    ginverse <- tryCatch(
      expr = {chol2inv(chol(G))},
      error = function(holder){return(NULL)}
    )

    if (is.null(ginverse)){
      stop("Matrix G is not positive definite.")
    }

    rownames(ginverse) <- rownames(G) # Add
    colnames(ginverse) <- colnames(G) # Add
    ginverse <- round(ginverse, digits)
  }
  if (rcn < rcn.thr & message){
    message('Reciprocal conditional number of G is: ',rcn, ', which is lower than the threshold of ',rcn.thr)
    message("G seems to be an ill-conditioned matrix.")
  }
  if (rcn < rcn.thr & !isTRUE(blend) & !isTRUE(bend) & !isTRUE(align)){
    stop('Consider bending, blending or aligning before invertion to make
         this matrix more stable.')
  }

  if (isTRUE(blend) || isTRUE(bend) || isTRUE(align)){
    Gb <- G.tuneup(G=G, A=A, blend=blend, pblend=pblend, bend=bend, align=align,
                   rcn=FALSE,sparseform=FALSE, determinant=FALSE, message=FALSE)$Gb

    # Try to get the inverse of Gb (tuned).
    ginverse <- tryCatch(
      expr = {chol2inv(chol(Gb))},
      error = function(holder){return(NULL)}
    )

    if (is.null(ginverse)){
      stop("Matrix G (tuned) is not positive definite.")
    }

    rownames(ginverse) <- rownames(Gb)
    colnames(ginverse) <- colnames(Gb)
  }

  ginverse <- round(ginverse, digits)

  # Reciprocal Condition Number RCN
  rcninv<-rcond(ginverse)
  if (message){
    message('Reciprocal conditional number for inverted matrix is: ', rcninv)
  }

  # Evaluating Matrix Empirically
  # Note: these are simple comparisons to check stability of matrix
  n <- ncol(ginverse)
  CN.1 <- ginverse[1,2]/sqrt(ginverse[1,1]*ginverse[1,1])
  CN.N <- ginverse[(n-1),n]/sqrt(ginverse[(n-1),(n-1)]*ginverse[n,n])
  max_diag <- abs(max(diag(ginverse)))
  max_off_diag <- max(abs(ginverse-diag(ginverse)))
  if(abs(CN.1)>0.99 | abs(CN.N)>0.99 | max_diag>1000 |  max_off_diag>1000) {
    if (message){
      message("Inverse of matrix G appears to be ill-conditioned.")
    }
    status <- 'ill-conditioned'
  } else {
    if (message){
      message("Inverse of matrix G does not appear to be ill-conditioned.")
    }
    status <- 'well-conditioned'
  }

  # Obtaining Matrix in Sparse Form if requested (ready for ASReml-R v.4)
  if (isTRUE(sparseform)) {
    ginverse.sparse <- full2sparse(ginverse)
    attr(ginverse.sparse, "INVERSE") <- TRUE
    return(list(Ginv.sparse = ginverse.sparse, rcn=rcninv, status=status))
  } else{
    attr(ginverse, "rowNames") <- rownames(ginverse)
    attr(ginverse, "colNames") <- colnames(ginverse)
    attr(ginverse, "INVERSE") <- TRUE
    return(list(Ginv = ginverse, rcn=rcninv, status=status))
  }

}
