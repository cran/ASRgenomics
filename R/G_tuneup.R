#' Tune-up the the genomic relationship matrix G
#'
#' @description
#' Generates a new matrix that can be \emph{blended}, \emph{bended} or \emph{aligned}
#' in order to make it stable for future use or matrix inversion. The input matrix should
#' be of the full form (\eqn{n \times n}) with individual names assigned to
#' \code{rownames} and \code{colnames}.
#'
#' This routine provides three options of tune-up:
#' \itemize{
#'  \item{\emph{Blend}. The \eqn{\boldsymbol{G}} matrix is blended (or averaged)
#'  with another matrix. \eqn{\boldsymbol{G}^\ast=(1-p) \boldsymbol{G} + p
#'  \boldsymbol{A}}, where \eqn{\boldsymbol{G}^\ast} is the blended matrix,
#'  \eqn{\boldsymbol{G}} is the original matrix,  \eqn{p} is the proportion of
#'  the identity matrix or pedigree-based relationship matrix to consider, and
#'  \eqn{\boldsymbol{A}} is the matrix to blend. Ideally, the pedigree-based
#'  relationship matrix should be used, but if this is not available (or it is
#'  of poor quality), then it is replaced by an identity matrix
#'  \eqn{\boldsymbol{I}}.}
#'  \item{\emph{Bend}. It consists on adjusting the original
#'  \eqn{\boldsymbol{G}} matrix to obtain a near positive definite matrix, which
#'  is done by making all negative or very small eigenvalues slightly positive.}
#'  \item{\emph{Align}. The original \eqn{\boldsymbol{G}} matrix is aligned to
#'  the matching pedigree relationship matrix \eqn{\boldsymbol{A}} where an
#'  \eqn{\alpha} and \eqn{\beta} parameters are obtained. More information can
#'  be found in the manual or in Christensen \emph{et al.} (2012).}
#' }
#'
#' The user should provide the matrices \eqn{\boldsymbol{G}} and
#' \eqn{\boldsymbol{A}} in full form (\eqn{n \times n})
#' and matching individual names should be
#' assigned to the \code{rownames} and \code{colnames} of the matrices.
#'
#' Based on procedures published by Nazarian and Gezan \emph{et al.} (2016).
#'
#' @param G Input of the genomic matrix \eqn{\boldsymbol{G}} to tune-up
#'  in full form (\eqn{n \times n}) (default = \code{NULL}).
#' @param A Input of the matching pedigree relationship matrix \eqn{\boldsymbol{A}}
#' in full form (\eqn{n \times n}) (default = \code{NULL}).
#' @param blend If \code{TRUE} a \emph{blending} with identity matrix \eqn{\boldsymbol{I}} or pedigree relationship matrix
#' \eqn{\boldsymbol{A}} (if provided) is performed (default = \code{FALSE}).
#' @param pblend If blending is requested this is the proportion of the identity matrix \eqn{\boldsymbol{I}} or
#' pedigree-based relationship matrix \eqn{\boldsymbol{A}} to blend for (default = \code{0.02}).
#' @param bend If \code{TRUE} a \emph{bending} is performed by making the matrix near positive definite (default = \code{FALSE}).
#' @param eig.tol Defines relative positiveness (\emph{i.e.}, non-zero) of eigenvalues compared to the largest one.
#' It determines which threshold of eigenvalues will be treated as zero (default = \code{1e-06}).
#' @param align If \code{TRUE} the genomic matrix \eqn{\boldsymbol{G}} is \emph{aligned} to the
#' matching pedigree relationship matrix \eqn{\boldsymbol{A}} (default = \code{FALSE}).
#' @param rcn If \code{TRUE} the reciprocal conditional number of the original and the bended,
#' blended or aligned matrix will be calculated (default = \code{TRUE}).
#' @param digits Set up the number of digits used to round the output matrix (default = \code{8}).
#' @param sparseform If \code{TRUE} it generates an inverse matrix in sparse form to be used directly in \pkg{asreml} with
#' required attributes (default = \code{FALSE}).
#' @param determinant If \code{TRUE} the determinant will be calculated, otherwise, this is obtained for
#' matrices of a dimension of less than 1,500 \eqn{\times} 1,500 (default = \code{TRUE}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with six of the following elements:
#' \itemize{
#' \item{\code{Gb}: the inverse of \eqn{\boldsymbol{G}} matrix in full form
#' (only if sparseform = \code{FALSE}).}
#' \item{\code{Gb.sparse}: if requested, the inverse of \eqn{\boldsymbol{G}} matrix in
#' sparse form (only if sparseform = \code{TRUE}).}
#' \item{\code{rcn0}: the reciprocal conditional number of the original matrix.
#' Values near zero are associated with an ill-conditioned matrix.}
#' \item{\code{rcnb}: the reciprocal conditional number of the blended, bended or aligned
#' matrix. Values near zero are associated with an ill-conditioned matrix.}
#' \item{\code{det0}: if requested, the determinant of the original matrix.}
#' \item{\code{blend}: if the matrix was \emph{blended}.}
#' \item{\code{bend}: if the matrix was \emph{bended}.}
#' \item{\code{align}: if the matrix was \emph{aligned}.}
#' }
#'
#' @references
#' Christensen, O.F., Madsen, P., Nielsen, B., Ostersen, T. and Su, G. 2012.
#' Single-step methods for genomic evaluation in pigs. Animal 6:1565-1571.
#' doi:10.1017/S1751731112000742.
#'
#' Nazarian A., Gezan S.A. 2016. GenoMatrix: A software package for pedigree-based
#' and genomic prediction analyses on complex traits. Journal of Heredity 107:372-379.
#'
#' @export
#'
#' @examples
#' # Example: Apple dataset.
#'
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#' G[1:5, 1:5]
#'
#' # Blend G matrix.
#' G_blended <- G.tuneup(G = G, blend = TRUE, pblend = 0.05)
#' G_blended$Gb[1:5, 1:5]
#'
#' # Bend G matrix.
#' G_bended <- G.tuneup(G = G, bend = TRUE, eig.tol = 1e-03)
#' G_bended$Gb[1:5, 1:5]
#'
#' \donttest{
#' # Example: Loblolly Pine dataset with pedigree - Aligned G matrix.
#'
#' A <- AGHmatrix::Amatrix(ped.pine)
#' dim(A)
#'
#' # Read and filter genotypic data.
#' M.clean <- qc.filtering(
#'  M = geno.pine655,
#'  maf = 0.05,
#'  marker.callrate = 0.2, ind.callrate = 0.20,
#'  na.string = "-9")$M.clean
#'
#' # Get G matrix.
#' G <- G.matrix(M = M.clean, method = "VanRaden", na.string = "-9")$G
#' G[1:5, 1:5]
#' dim(G)
#'
#' # Match G and A.
#' Aclean <- match.G2A(A = A, G = G, clean = TRUE, ord = TRUE, mism = TRUE)$Aclean
#'
#' # Align G with A.
#' G_align <- G.tuneup(G = G, A = Aclean, align = TRUE)
#' G_align$Gb[1:5, 1:5]
#' }
#'

G.tuneup <- function(G = NULL, A = NULL, blend = FALSE, pblend = 0.02,
                     bend = FALSE, eig.tol = 1e-06, align = FALSE, rcn = TRUE,
                     digits = 8, sparseform = FALSE, determinant = TRUE, message = TRUE){

  # Check if the class of G is matrix
  if (is.null(G) || !inherits(G, "matrix")) {
    stop('G should be a valid object of class matrix.')
  }
  # Check the rownames/colnames of G
  if (is.null(rownames(G))){
    stop('Individual names not assigned to rows of matrix G.')
  }
  if (is.null(colnames(G))){
    stop('Individual names not assigned to columns of matrix G.')
  }
  if ((identical(rownames(G), colnames(G))) == FALSE){
    stop("Rownames and colnames of matrix G do not match.")
  }
  # Checks on other parameters
  if (pblend < 0 | pblend > 1) {
    stop("Specification of pblend must be between 0 and 1.")
  }
  if (eig.tol <= 0) {
    stop("Value for eig.tol must be positive.")
  }

  # Check option parameters
  if (isTRUE(blend) & isTRUE(bend) & isTRUE(align)){
    stop('More than one option (blend, bend or align) was requested. Choose only one of them')
  }
  if (!isTRUE(blend) & !isTRUE(bend) & !isTRUE(align)){
    stop('No option (blend, bend or aling) were requested.')
  }
  n <- dim(G)[1]
  p <- dim(G)[2]

  # #G <- as.matrix(forceSymmetric(G))    #Force a square matrix to be a symmetric Matrix
  # if (p != n) { stop('Matrix G is not a square matrix.') }

  if (!is.null(A)) {
    # Check if the class of A is matrix
    if (!inherits(A, "matrix")) {
      stop('A should be a valid object of class matrix.')
    }
    if (is.null(rownames(A))){
      stop('Individual names not assigned to rows of matrix A.')
    }
    if (is.null(colnames(A))){
      stop('Individual names not assigned to columns of matrix A.')
    }
    nA <- dim(A)[1]
    pA <- dim(A)[2]
    #A <- as.matrix(forceSymmetric(A))    #Force a square matrix to be a symmetric Matrix
    if (pA != nA) { stop('Matrix A is not a square matrix.') }
    if (p != pA | n != nA) { warning('Matrix A and G are not of the same dimensions.
                                     Use match_G2A() to obtain subset of matched genotypes.') }
    }

  # Reciprocal Condition Number RCN
  if (isTRUE(rcn)){
  rcn0 <- rcond(G)
  if (message) {
    message('Reciprocal conditional number for original matrix is: ', rcn0)
  }
  } else { rcn0 <- NULL}

  # Checking for determinant
  if (n > 1500) {
    if (!determinant) {
      if (message) {
        message('Determinant is not calculated as matrix is too large.')
      }
      detG <- NULL
    } else {
      detG <- det(G)
      if (message) {
        message('Determinant for original matrix is: ', detG)
      }
    }
  } else {
    detG <- det(G)
    if (message) {
      message('Determinant for original matrix is: ', detG)
    }
  }



  # Performing blend with I if requested
  if (isTRUE(blend)){
    if (is.null(A)) {
      Gb <- (1-pblend)*G + pblend*diag(x=1, nrow=n, ncol=n)
      if (message) {
        message('Matrix was BLENDED using an identity matrix.')
      }
    }
    # Performing blend with A if requested
    if (!is.null(A)) {
      if (p != pA | n != nA) { stop('Matrix A and G are not of the same dimension.') }
      if ( sum(rownames(G) == rownames(A)) < n) {
        warning('Names of rows/columns do not match between G and A matrix.')
        stop('You should use the function match.G2A to match these two matrices.')
      }
      if ( sum(rownames(G) == rownames(A)) == n) {
        Gb <- (1-pblend)*G + pblend*A
        if (message) {
          message('Matrix was BLENDED using the provided A matrix.')
        }
      }
    }
  }

  # Obtaining the near positive definite matrix (bend)
  if (isTRUE(bend)) {
    Gb <- as.matrix(Matrix::nearPD(G, corr=FALSE, keepDiag=FALSE, do2eigen=TRUE,
                           doSym=TRUE, doDykstra=TRUE, only.values=FALSE,
                           eig.tol=eig.tol, conv.tol=1e-07, posd.tol=1e-02,
                           maxit=100, conv.norm.type="I", trace=FALSE)$mat)
    if (message) {
      message('Matrix was BENDED.')
    }
  }

  # Performing alignment with A matrix (align)
  if (isTRUE(align)) {
    if (is.null(A)) {
    stop('A matrix needs to be provided for align be performed.')}
    Xm <- cbind(c(1,1),c(mean(diag(G)),mean(G)))
    y <- as.matrix(c(mean(diag(A)),mean(A)))
    beta <- solve(Xm) %*% y
    Gb <- beta[1] + beta[2]*G
    if (message) {
      message('Matrix was ALIGNED.')
    }
  }

  # Rounding Gb
  Gb <- round(Gb, digits)
  ### Need to make sure rownames(Gb) and colnames(Gb) are correct

  # Reciprocal Condition Number RCN
  if (isTRUE(rcn)){
  rcnb <- rcond(Gb)
  if (message) {
    message('Reciprocal conditional number for tune-up matrix is: ', rcnb)
  }
  } else { rcnb <- NULL}

  # Obtaining Matrix in Sparse Form if requested (ready for ASReml-R v.4)
  if (isTRUE(sparseform)) {
    Gb.sparse <- full2sparse(Gb)
    return(list(Gb.sparse=Gb.sparse, rcn0=rcn, det0=detG, rcnb=rcnb, blend=blend, bend=bend, align=align))
  } else {
    return(list(Gb=Gb, rcn0=rcn0, det0=detG, rcnb=rcnb, blend=blend, bend=bend, align=align))
  }


}
