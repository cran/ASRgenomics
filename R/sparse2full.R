#' Generates a full matrix form from a sparse form matrix
#'
#' Modifies the input sparse form matrix into its full form.
#' The sparse form has three columns per line, corresponding to the set:
#' \code{Row, Col, Value}, and is defined by a lower triangle row-wise
#' of the full matrix and is sorted as columns within row.
#' Individual names should be assigned as attributes: \code{attr(K, "rowNames")}
#' and \code{attr(K, "colNames")}. If these are not provided they are considered
#' as 1 to \eqn{n}.
#'
#' Based on a function from ASReml-R 3 library by Butler \emph{et al.} (2009).
#'
#' @param K A square matrix in sparse form (default = \code{NULL}).
#'
#' @return A full square matrix where individual names are assigned to
#' \code{rownames} and \code{colnames}.
#' If attribute \code{INVERSE} is found this is also passed to the full matrix.
#'
#' @references
#' Butler, D.G., Cullis, B.R., Gilmour, A.R., and Gogel, B.J. 2009.
#' ASReml-R reference manual. Version 3. The Department of Primary
#' Industries and Fisheries (DPI&F).
#'
#' @export
#'
#' @examples
#' # Get G matrix.
#' Gsp <- G.matrix(M = geno.apple, method = "VanRaden", sparseform = TRUE)$G.sparse
#' head(Gsp)
#' head(attr(Gsp, "rowNames"))
#'
#' # Transform into full matrix.
#' G <- sparse2full(K = Gsp)
#' G[1:5, 1:5]
#'

sparse2full <- function(K = NULL) {

  if (is.null(K) || !inherits(K, "matrix")) {
    stop('K should be a valid object of class matrix.')
  }

  # Collect inverse attribute if any.
  INVERSE <- attr(K, "INVERSE")

  # Collect dimnames to apply on new matrix.
  rownames <- attr(K, "rowNames")

  # Collect number of rows.
  nrow <- max(K[, 1])

  # Collect number of columns.
  ncol <- max(K[, 2])

  # Get dummy rownames if not provided in attributes.
  if (is.null(rownames)) {
    rownames <- as.character(1:nrow)
  }

  # Assign relationships to relative positions.
  K.full <- rep(0, nrow * ncol)
  K.full[(K[, 2] - 1) * nrow + K[, 1]] <- K[, 3]
  K.full[(K[, 1] - 1) * nrow + K[, 2]] <- K[, 3]

  # Reshape to matrix.
  K.full <- matrix(K.full, nrow = nrow, ncol = ncol, byrow = FALSE)

  # Assign row and colnames.
  attr(K.full, "colNames") <-
    attr(K.full, "rowNames") <-
      colnames(K.full) <-
        rownames(K.full) <-
          rownames

  if (!is.null(INVERSE)) {attr(K.full, "INVERSE") <- INVERSE}

  return(K.full)

}
