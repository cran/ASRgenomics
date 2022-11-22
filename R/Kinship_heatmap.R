#' Enhanced heatmap plot for a kinship matrix K
#'
#' Generates a heatmap with dendrogram based on a provided kinship matrix.
#' This matrix can be a pedigree relationship matrix \eqn{\boldsymbol{A}}, a
#' genomic relationship matrix \eqn{\boldsymbol{G}} or a hybrid relationship
#' matrix \eqn{\boldsymbol{H}}.
#' Individual names should be assigned to \code{rownames} and \code{colnames}.
#' It sorts individuals according to dendrogram in both columns and rows.
#'
#' Uses the library \code{superheat} from Barter and Yu (2018) to generate plots.
#'
#' @param K Input of a kinship matrix in full format (\eqn{n \times n}) (default = \code{NULL}).
#' @param dendrogram If \code{TRUE} a dendrogram is added to the columns based on the
#' kinship matrix (default = \code{TRUE}).
#' @param clustering.method The clustering method considered for the dendrogram.
#' Options are: \code{"hierarchical"} and \code{"kmeans"} (default = \code{"hierarchical"}).
#' @param dist.method The method considered to calculate the distance matrix between
#' individuals used for hierarchical clustering. Options are: \code{"euclidean"},
#' \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"} and
#' \code{"minkowski"} (default = \code{"euclidean"}).
#' @param row.label If \code{TRUE} the individual names (\code{rownames}) are added as labels to
#' the left of the heatmap (default = \code{TRUE}).
#' @param col.label If \code{TRUE} the individual names (\code{colnames}) are added as labels to
#' the bottom of the heatmap (default = \code{FALSE}).
#'
#' @return
#' A plot with the properties specified by the above arguments.
#'
#' @references
#' Barter, R.L. and Yu, B. 2018. Superheat: An R package for creating beautiful
#' and extendable heatmaps for visualizing complex data.
#' J. Comput. Graph. Stat. 27(4):910-922.
#'
#' @export
#'
#' @examples
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#' G[1:5, 1:5]
#'
#' # Plot a subset of the individuals.
#' kinship.heatmap(K = G[1:10, 1:10], dendrogram = TRUE, row.label = TRUE, col.label = TRUE)
#'


kinship.heatmap <- function(K = NULL, dendrogram = TRUE,
                            clustering.method = c("hierarchical", "kmeans"),
                            dist.method = c("euclidean", "maximum", "manhattan",
                                            "canberra", "binary", "minkowski"),
                            row.label = TRUE, col.label = FALSE){

  # Check if the class of K is matrix
  if (is.null(K) || !inherits(K, "matrix")) {
    stop("K should be a valid object of class matrix.")
  }
  # Check the rownames/colnames
  if (is.null(rownames(K))){
    stop("Individual names not assigned to rows of matrix K.")
  }
  if (is.null(colnames(K))){
    stop('Individual names not assigned to columns of matrix K.')
  }
  if ((identical(rownames(K), colnames(K))) == FALSE){
    stop("Rownames and colnames of matrix K do not match.")
  }
  # Check if the are missing values
  if (any(is.na(K))){
    stop("Matrix K contains some missing data.")
  }
  clustering.method <- match.arg(clustering.method)
  dist.method <- match.arg(dist.method)

  if (isTRUE(col.label)) {
    labCol <- "variable"
  } else {
    labCol <- "none"
  }
  if (isTRUE(row.label)) {
    labRow <- "variable"
  } else {
    labRow <- "none"
  }

  pp <- superheat(K,
            pretty.order.rows = TRUE,
            pretty.order.cols = TRUE,
            col.dendrogram = dendrogram,
            dist.method = dist.method,
            clustering.method = clustering.method,
            scale = FALSE,
            left.label = labRow,
            left.label.col = "white",
            force.left.label = TRUE,
            bottom.label = labCol,
            bottom.label.col = "white",
            bottom.label.text.angle = 90,
            bottom.label.text.size = 2.5,
            force.bottom.label = TRUE,
            print.plot = TRUE,
            legend.text.size = 10, # default 12
            legend.width = 2.0,
            left.label.text.size = 2.5)

  return(pp)
}


