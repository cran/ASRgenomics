#' Generates the inverse of the hybrid H matrix
#'
#' The single-step GBLUP approach combines the information from the pedigree relationship matrix
#' \eqn{\boldsymbol{A}} and the genomic relationship matrix \eqn{\boldsymbol{G}} in one
#' hybrid relationship matrix called \eqn{\boldsymbol{H}}.
#' This function will calculate directly the inverse of this matrix \eqn{\boldsymbol{H}}.
#' The user should provide the matrices \eqn{\boldsymbol{A}} or
#' its inverse (only one of these is required) and the
#' inverse of the matrix \eqn{\boldsymbol{G}} (\eqn{\boldsymbol{G_{inv}}}) in its full form. Individual names should
#' be assigned to  \code{rownames} and \code{colnames}, and individuals from
#' \eqn{\boldsymbol{G_{inv}}} are verified to be all a subset within individuals from
#' \eqn{\boldsymbol{A}} (or \eqn{\boldsymbol{A_{inv}}}).
#'
#' @param A Input of the pedigree relationship matrix \eqn{\boldsymbol{A}}
#' in full form (\eqn{na \times na}) (default = \code{NULL}).
#' @param Ainv Input of the inverse of the pedigree relationship matrix
#' \eqn{\boldsymbol{A}^{-1}} in full form (\eqn{na \times na}) (default = \code{NULL}).
#' @param Ginv Input of the inverse of the genomic relationship matrix
#' \eqn{\boldsymbol{G}^{-1}} in full form (\eqn{ng \times ng}) (default = \code{NULL}).
#' @param lambda The scaling factor for \eqn{(\boldsymbol{G}^{-1}-\boldsymbol{A}^{-1}_{22})} (default = \code{NULL}).
#' @param tau The scaling factor for \eqn{\boldsymbol{G}^{-1}} (default = \code{1}).
#' @param omega The scaling factor for \eqn{\boldsymbol{A}^{-1}_{22}} (default = \code{1}).
#' @param sparseform If \code{TRUE} it generates the requested matrix in sparse form to be used
#' directly in \pkg{asreml} with required attributes (default = \code{FALSE}).
#' @param keep.order If \code{TRUE} the original order of the individuals from the
#' \eqn{\boldsymbol{A}} or \eqn{\boldsymbol{A_{inv}}} matrix is kept.
#' Otherwise the non-genotyped individuals are placed first and
#' then genotyped individuals (default = \code{TRUE}).
#' @param digits Set up the number of digits used to round the output matrix (default = \code{8}).
#' @param inverse If \code{TRUE} it generates the inverse of \eqn{\boldsymbol{H}} matrix (default = \code{TRUE})
#' (to be deprecated).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return
#' The inverse of the hybrid matrix \eqn{\boldsymbol{H}} matrix, in full or sparse form with
#' required attributes to be used in \pkg{asreml}.
#'
#' @details
#' The generation of the \eqn{\boldsymbol{H^{-1}}} matrix contains a few scaling factors
#' to help with the calculation of this inverse and
#' to allow further exploration of the combination of the
#' information from the \eqn{\boldsymbol{A^{-1}}} and \eqn{\boldsymbol{G^{-1}}}.
#' We follow the specifications described by Martini \emph{et. al} (2018),
#' which is done by specifying the parameters \eqn{\lambda}, or the pair
#' \eqn{\tau} and \eqn{\omega}.
#'
#' \if{html}{
#' The general expression used is:
#'   \deqn{\boldsymbol{H^{-1}}=\boldsymbol{A^{-1}}+\begin{bmatrix}\boldsymbol{0}&\boldsymbol{0}\\\boldsymbol{0}&(\tau\boldsymbol{G^{-1}}-\omega\boldsymbol{{A_{22}^{-1}}})\end{bmatrix}}
#' }
#'
#' \if{html}{
#' and a more common representation of the above expression is found when \eqn{\tau = \omega = \lambda}, as shown below:
#'   \deqn{\boldsymbol{H^{-1}}=\boldsymbol{A^{-1}}+\begin{bmatrix}\boldsymbol{0}&\boldsymbol{0}\\\boldsymbol{0}&\lambda(\boldsymbol{G^{-1}}-\boldsymbol{{A_{22}^{-1}}})\end{bmatrix}}
#' }
#'
#' If \code{inverse = FALSE} the \eqn{\boldsymbol{H}}
#' matrix is provided instead of its inverse. This option will be deprecated and
#' it is better to use the function \link{H.matrix}.
#'
#' \if{html}{
#' The \eqn{\boldsymbol{H}} matrix is obtained with the following equations:
#'   \deqn{\boldsymbol{H}=\boldsymbol{A}+\begin{bmatrix}\boldsymbol{A}_{12}\boldsymbol{A}_{22}^{-1}(\boldsymbol{G}-\boldsymbol{A}_{22})\boldsymbol{A}_{22}^{-1}\boldsymbol{A}_{21}&\boldsymbol{A}_{12}\boldsymbol{A}_{22}^{-1}(\boldsymbol{G}-\boldsymbol{A}_{22})\\(\boldsymbol{G}-\boldsymbol{A}_{22})\boldsymbol{A}_{22}^{-1}\boldsymbol{A}_{21}&(\boldsymbol{G}-\boldsymbol{A}_{22})\end{bmatrix}}
#' }
#'
#' @references
#' Christensen, O.F., Lund, M.S. 2010. Genomic prediction matrix when some animals
#' are not genotyped. Gen. Sel. Evol. 42(2):1–8.
#'
#' Christensen, O., Madsen, P., Nielsen, B., Ostersen, T., and Su, G. 2012. Single-step methods
#' for genomic evaluation in pigs. Animal 6(10):1565–1571.
#'
#' Legarra, A., Aguilar, I., and Misztal, I. 2009. A relationship matrix including full
#' pedigree and genomic information. J. Dairy Sci. 92:4656-4663.
#'
#' Martini, J.W.R., Schrauf, M.F., Garcia-Baccino, C.A., Pimentel, E.C.G., Munilla, S.,
#' Rogberg-Muñoz, A., Cantet, R.J.C., Reimer, C., Gao, N., Wimmer, V., and Simianer, H. 2018.
#' The effect of the \eqn{H^{-1}} scaling factors \eqn{\tau} and \eqn{\omega}
#' on the structure of \eqn{H} in the single-step procedure.
#' Genet. Sel. Evol. 50:1-9.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Get A matrix.
#' A <- AGHmatrix::Amatrix(data = ped.pine)
#' A[1:5,1:5]
#' dim(A)
#'
#' # Read and filter genotypic data.
#' M.clean <- qc.filtering(
#'  M = geno.pine655,
#'  maf = 0.05,
#'  marker.callrate = 0.2, ind.callrate = 0.20,
#'  na.string = "-9",
#'  plots = FALSE)$M.clean
#'
#' # Get G matrix.
#' G <- G.matrix(M.clean, method = "VanRaden", na.string = "-9")$G
#' G[1:5, 1:5]
#' dim(G)
#'
#' # Match G and A.
#' check <- match.G2A(
#'  A = A, G = G,
#'  clean = TRUE, ord = TRUE, mism = TRUE, RMdiff = TRUE)
#'
#' # Align G matrix with A.
#' G_align <- G.tuneup(G = check$Gclean, A = check$Aclean, align = TRUE, sparseform = FALSE)$Gb
#'
#' # Get Ginverse using the G aligned.
#' Ginv <- G.inverse(G = G_align, sparseform = FALSE)$Ginv
#' Ginv[1:5, 1:5]
#' dim(Ginv)
#'
#' # Obtain Hinv.
#' Hinv <- H.inverse(A = A, Ginv = Ginv, lambda = 0.90, sparseform = TRUE)
#' head(Hinv)
#' attr(Hinv, "INVERSE")
#' }
#'

H.inverse <- function(A = NULL, Ainv = NULL, Ginv = NULL,
                      lambda = NULL,tau = 1, omega = 1,
                      sparseform = FALSE, keep.order = TRUE,
                      digits = 8, inverse = TRUE,
                      message = TRUE){

  # Check if we have lambda and make some checks
  if(!is.null(lambda)){
    tau <- lambda
    omega <- lambda
    if(message){
      message("A lambda value was provided and it will be used instead of tau and omega.")
    }
    if (lambda < 0 | lambda > 1) {
      stop("Value of lambda must be between 0 and 1.")
    }
  }else{
    if(message){
      message("No lambda value was provided, tau and omega scaling factors will be considered.")
    }
  }


  # Checks on some input parameters
  if (tau < 0 | tau > 2) {
    stop("Value of tau must be between 0 and 1.")
  }
  if (omega < -1 | omega > 1) {
    stop("Value of omega must be between -1 and 1.")
  }
  if (tau == 0 & omega == 0) {
    stop("Values of tau and omega can not be both equal to zero.")
  }


  # We need A or Ainv
  if (is.null(A) & is.null(Ainv)) {
    stop('A or Ainv needs to be provided.')
  }

  if (!is.null(Ainv)) { # We have Ainv and get A
    if (!inherits(Ainv, "matrix")) {
      stop("Ainv should be a valid object of class matrix.")
    }
    if (is.null(rownames(Ainv))){
      stop("Individual names not assigned to rows of matrix Ainv.")
    }
    if (is.null(colnames(Ainv))){
      stop("Individual names not assigned to columns of matrix Ainv.")
    }
    if ((identical(rownames(Ainv), colnames(Ainv))) == FALSE){
      stop("Rownames and colnames of matrix Ainv do not match.")
    }
    crit<-sum(diag(Ainv)>=2)
    if (crit == 0) {
      warning('It seems like an A matrix was provided not an Ainv matrix.')
    }
    A <- G.inverse(G=Ainv, blend=FALSE, bend=FALSE, message = FALSE)$Ginv
    attributes(A)$INVERSE <- NULL
  } else { # Then we have A and we get Ainv to use
    if (!inherits(A, "matrix")) {
      stop("A should be a valid object of class matrix.")
    }
    if (is.null(rownames(A))){
      stop("Individual names not assigned to rows of matrix A.")
    }
    if (is.null(colnames(A))){
      stop("Individual names not assigned to columns of matrix A.")
    }
    if ((identical(rownames(A), colnames(A))) == FALSE){
      stop("Rownames and colnames of matrix A do not match.")
    }
    crit<-sum(diag(A)>2)
    if (crit >0) {
      warning('It seems like an Ainv matrix was provided not an A matrix.')
    }
    Ainv <- G.inverse(G=A, blend=FALSE, bend=FALSE, message = FALSE)$Ginv

  }


  # We need Ginv
  # Check components of Ginv
  if (is.null(Ginv) || !inherits(Ginv, "matrix")) {
    stop("Ginv should be a valid object of class matrix.")
  }
  if (is.null(rownames(Ginv))){
    stop("Individual names not assigned to rows of matrix Ginv.")
  }
  if (is.null(colnames(Ginv))){
    stop("Individual names not assigned to columns of matrix Ginv.")
  }
  if ((identical(rownames(Ginv), colnames(Ginv))) == FALSE){
    stop("Rownames and colnames of matrix Ginv do not match.")
  }

  # Check for consistency between Ainv and Ginv
  Aind <- rownames(Ainv)
  Gind <- rownames(Ginv)

  notGenotyped <- which((Aind %in% Gind) == FALSE) # which ind has Ainv but not Ginv
  notPedigree <- which(Gind %in% Aind == FALSE)    # which ind has Ginv but no Ainv

  # If Ginv has different ind than Ainv, stop here.
  if (length(notPedigree) > 0) {
    stop("Matrix Ginv has ", length(notGenotyped), " individuals that are not present in A or Ainv.")
  }
  # Check if Ainv has different ind than Ginv
  if (length(notGenotyped) > 0 & message){
    message("Matrix A or Ainv has ", length(notGenotyped), " individuals that are not present in Ginv.")
  }
  if (length(notGenotyped) < 0 & message){
    message("Matrix A or Ainv has all individuals that are present in Ginv.")
  }


  # Creating indexes (sorted matrices)
  idA <- rownames(Ainv)
  idG <- rownames(Ginv)
  idH <- unique(c(idG,idA))
  idH <- rev(idH)
  A <- A[idH,idH]
  Ainv <- Ainv[idH,idH]
  genotyped <- idH %in% idG == TRUE
  A11inv <- Ainv[!genotyped, !genotyped]
  A12inv <- Ainv[!genotyped, genotyped]
  A22inv <- solve(A[genotyped,genotyped])
  Ginv <- Ginv[idH[genotyped],idH[genotyped]]

  # Making sure both Ginv and A22inv agree on sorted individuals
  if (all(rownames(Ginv) == row.names(A22inv)) == FALSE) {  # Then we have a problem
    stop('Order of matrix Ginv does not match subset of A or Ainv.')
  }

  # Traditional method using Ginv & Ainv to obtain Hinv directly.
  if (isTRUE(inverse)){
    # Obtaining Hinv
    H11inv <- matrix(0, nrow=dim(A11inv)[1], ncol=dim(A11inv)[2] )
    H12inv <- matrix(0, nrow=dim(A12inv)[1], ncol=dim(A12inv)[2] )
    H22inv <- tau*Ginv - omega*A22inv
    Hinv <- Ainv + cbind(rbind(H11inv, t(H12inv)), rbind(H12inv, H22inv))

    # round Hinv
    Hinv <- round(Hinv, digits)

    # Same order as Ainverse (if requested)
    if (keep.order) {
      Hinv <- Hinv[match(Aind, rownames(Hinv)), match(Aind, rownames(Hinv))]
    }
    attr(Hinv, "INVERSE") <- TRUE

    # Generating sparseform matrix
    if (sparseform) {
      Hinv <- full2sparse(Hinv)
      attr(Hinv, "INVERSE") <- TRUE
    }
    return(Hinv=Hinv)
  } else {
    A11 <- A[!genotyped, !genotyped]
    A12 <- A[!genotyped, genotyped]
    A21 <- t(A12)
    A22 <- A[genotyped,genotyped]

    #H22 = solve((tau * Ginv + (1 - omega) * A22inv))
    H22 <- G.inverse((tau * Ginv + (1 - omega) * A22inv), message = FALSE)
    if (H22$status == 'ill-conditioned') {
      stop("Matrix H is ill-conditioned, try different scaling factors.")
    } else {
      H22  <- H22$Ginv
    }

    # Old repetitive code.
    # H11 = A12 %*% A22inv %*% (H22 - A22) %*% A22inv %*% A21
    # H12 = A12 %*% A22inv %*% (H22 - A22)
    # H21 = (H22 - A22) %*% A22inv %*% A21
    # H22 = (H22 - A22)

    # Pre-calculated matrices.
    A22A21 <- A22inv %*% A21

    # Finalize multiplications.
    H22 <- (H22 - A22)
    H12 <- A12 %*% A22inv %*% H22
    H11 <- H12 %*% A22A21
    H21 <- H22 %*% A22A21

    # Bind matrices.
    H <- A + cbind(rbind(H11, H21), rbind(H12, H22))

    # round H
    H <- round(H, digits)

    # Same order as A (if requested)
    if (keep.order) {
      H <- H[match(Aind, rownames(H)), match(Aind, rownames(H))]
    }

    # Generating sparseform matrix
    if (sparseform) {
      H <- full2sparse(H)
      attr(H, "INVERSE") <- FALSE
    }
    return(H=H)
  }

}
