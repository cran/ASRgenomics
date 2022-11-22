#' Generates the conditional predictions of random effects (BLUPs)
#'
#' @description
#' Predicts random effects values for individuals with unobserved responses (here called \code{x},
#' a vector of length \eqn{nx}) based on known random effect values for individuals with
#' observed responses (here called \code{y}, a vector of length \eqn{ny}). This is done using the
#' common genomic relationship matrix  \eqn{\boldsymbol{G}} for all
#' individuals (full matrix of dimension \eqn{(nx + ny) \times (nx + ny)}).
#'
#' The prediction of unobserved responses will be performed through the
#' multivariante Normal conditional distribution. These predictions are identical to
#' what would be obtained if the entire set of individuals (\eqn{nx + ny}) were included into a
#' GBLUP animal model fit with individuals in the set \code{x} coded as missing.
#'
#' The user needs to provide the matrix \eqn{\boldsymbol{G}} in full form.
#' Individual names (\eqn{nx + ny}) should be assigned to \code{rownames} and \code{colnames}, and these
#' can be in any order. If the variance-covariance matrix of the set \code{y} is provided,
#' standard errors of random effects in set \code{x} are calculated.
#'
#' @param G Input of the genomic relationship matrix \eqn{\boldsymbol{G}} in full form
#' (of dimension \eqn{(nx + ny) \times (nx + ny)}) (default = \code{NULL}).
#' @param gy Input of random effects (\emph{e.g.} breeding values) for individuals with known values.
#' Individual names should be assigned to \code{rownames} of this vector
#' and be found on the matrix \eqn{\boldsymbol{G}} (default = \code{NULL}).
#' @param vcov.gy The variance-covariance matrix associated with the random effects from the
#' individuals with known values (set \code{y}, of dimension \eqn{ny \times ny})
#' (default = \code{NULL}).
#'
#' @return A data frame with the predicted random effect values for individuals with
#' unobserved responses in the set \code{x}. If the variance-covariance matrix is provided,
#' standard errors are included in an additional column.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(asreml) # Load asreml.
#'
#' # Example: Apple data creating 100 missing observations.
#'
#' # Prepare G (nx + ny).
#' G <- G.matrix(M = geno.apple, method = "VanRaden", sparseform = FALSE)$G
#' dim(G)
#'
#' # Prepare subset of data.
#' # Select only 147 from 247 individuals from pheno.apple and geno.apple.
#' Gy <- G[1:147, 1:147]
#' phenoy <- pheno.apple[1:147, ]
#'
#' # Obtain the BLUPs for the 147 individuals using ASReml-R.
#'
#' # Blend Gy.
#' Gyb <- G.tuneup(G = Gy, blend = TRUE, pblend = 0.02)$Gb
#'
#' # Get the Gy inverse
#' Gyinv  <- G.inverse(G = Gyb, sparseform = TRUE)$Ginv.sparse
#'
#' # Fit a GBLUP model
#' phenoy$INDIV <- as.factor(phenoy$INDIV)
#' modelGBLUP <-
#'  asreml(
#'   fixed = JUI_MOT ~ 1,
#'   random = ~vm(INDIV, Gyinv),
#'   workspace = 128e06,
#'   data = phenoy)
#'
#' # Obtain Predictions - BLUP (set y).
#' BLUP <- summary(modelGBLUP,coef = TRUE)$coef.random
#' head(BLUP)
#' gy <- as.matrix(BLUP[, 1])
#' rownames(gy) <- phenoy$INDIV
#'
#' # Ready to make conditional predictions.
#' g.cond <- G.predict(G = G, gy = gy)
#' head(g.cond)
#' }
#'

G.predict <- function(G = NULL, gy = NULL, vcov.gy = NULL){

  # Check if the class of G is matrix
  if (is.null(G) || !inherits(G, "matrix")) {
    stop("G should be a valid object of class matrix.")
  }
  if (is.null(rownames(G))){
    stop("Individual names not assigned to rows of matrix G.")
  }
  if (is.null(colnames(G))){
    stop("Individual names not assigned to columns of matrix G.")
  }
  if ((identical(rownames(G), colnames(G))) == FALSE){
    stop("Rownames and colnames of matrix G do not match.")
  }

  # Check if the class of gy is matrix
  if (is.null(gy) || !inherits(gy, "matrix")) {
    stop("gy should be a valid object of class matrix.")
  }
  if (is.null(rownames(gy))){
    stop("Individual names not assigned to rows of gy.")
  }

  # Check for consistency between G and resp
  resp <- gy
  respind <- rownames(resp)
  Gind <- rownames(G)

  yidx <-  which(Gind %in% respind == TRUE)
  xidx <- which(Gind %in% respind == FALSE)
  if (length(xidx) == 0) {
    stop("Individuals in G are the same as those in the vector gy, nothing to predict.")
  }

  # Creating indexes
  G_yy <- G[yidx, yidx]
  G_xy <- G[xidx, yidx]

  # Obtain inverse of G_yy
  G_yy.inv <- G.inverse(G=G_yy, bend=FALSE, blend=FALSE, sparseform=FALSE,
                        message=FALSE)
  if (G_yy.inv$status == 'ill-conditioned') {
    stop("Inverse of matrix G is ill-conditioned, use G.tuneup and provide G.")
  } else {
    G_yy.inv  <- G_yy.inv$Ginv
  }

  # Obtain predictions
  beta <- G_xy %*% G_yy.inv
  pred <- beta %*% resp
  pred <- data.frame(predicted.value=pred)

  # Obtaining vcov(preds)
  if (!is.null(vcov.gy)) {
    # Check if the class of Kinv is matrix
    if (is.null(vcov.gy) || !inherits(vcov.gy, "matrix")) {
      stop("Input vcov.resp should be a valid object of class matrix.")
    }
    if (is.null(rownames(vcov.gy))){
      stop("Individual names not assigned to rows of matrix vcov.gy.")
    }
    if (is.null(colnames(vcov.gy))){
      stop("Individual names not assigned to columns of matrix vcov.gy.")
    }
    if ((identical(row.names(vcov.gy), colnames(vcov.gy))) == FALSE){
      stop("Rownames and colnames of matrix vcov.gy do not match.")
    }
    vcov.pred <- beta %*% vcov.gy %*% t(beta)
    st.pred <- sqrt(diag(vcov.pred))
    pred <- data.frame(predicted.value=pred, std.error=st.pred)
  }

  return(pred=pred)
}
