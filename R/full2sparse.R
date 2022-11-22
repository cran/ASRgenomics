#' Generates a sparse form matrix from a full form matrix
#'
#' Modifies the input square matrix into its sparse form with three columns per line,
#' corresponding to the set: \code{Row, Col, Value}. This matrix defines the lower triangle
#' row-wise of the original matrix and it is sorted as columns within row.
#' Values of zero on the matrix are dropped by default. Individual
#' names should be assigned to \code{rownames} and \code{colnames}.
#'
#' Based on the function published by Borgognone \emph{et al.} (2016).
#'
#' @param K A square matrix in full form (\eqn{n \times n}) (default = \code{NULL}).
#' @param drop.zero If \code{TRUE} observations equal to zero are dropped
#' (default = \code{TRUE}).
#'
#' @return A matrix in sparse form with columns: \code{Row, Col, Value}
#' together with the attributes \code{rowNames} and \code{colNames}.
#' If attribute \code{INVERSE} is found this is also passed to the sparse matrix.
#'
#' @references
#' Borgognone, M.G., Butler, D.G., Ogbonnaya, F.C and Dreccer, M.F. 2016.
#' Molecular marker information in the analysis of multi-environment trial helps
#' differentiate superior genotypes from promising parents. Crop Science 56:2612-2628.
#'
#' @export
#'
#' @examples
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden", sparseform = FALSE)$G
#' G[1:5, 1:5]
#'
#' # Transform matrix into sparse.
#' G.sparse <- full2sparse(K = G)
#' head(G.sparse)
#' head(attr(G.sparse, "rowNames"))
#' head(attr(G.sparse, "colNames"))
#'


full2sparse <- function (K = NULL, drop.zero = TRUE) {

  if (is.null(K) || !inherits(K, "matrix")) {
    stop('K should be a valid object of class matrix.')
  }
  # Check the attributes of K
  if (is.null(rownames(K))){
    stop('Individual names not assigned to rows of matrix K.')
  }
  if (is.null(colnames(K))){
    stop('Individual names not assigned to columns of matrix K.')
  }
  if ((identical(rownames(K), colnames(K))) == FALSE){
    stop("Rownames and colnames of matrix K do not match.")
  }

  INVERSE <- attr(K, "INVERSE")

  if(drop.zero) {
    which <- (K != 0 & lower.tri(K, diag = TRUE))
  } else {
    which <- lower.tri(K, diag = TRUE)
  }

  sparse.frame <- data.frame(Row = t(row(K))[t(which)], Col = t(col(K))[t(which)],
                   Value = t(K)[t(which)])

  sparse.frame <- as.matrix(sparse.frame)

  # Add attributes.
  attr(sparse.frame, "rowNames") <- rownames(K)
  attr(sparse.frame, "colNames") <- colnames(K)

  # Pass collected attributes.
  if (!is.null(INVERSE)) {attr(sparse.frame, "INVERSE") <- INVERSE}

  return(sparse.frame)

}
