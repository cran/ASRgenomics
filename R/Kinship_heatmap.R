#' Build dendrogram segment data for ggplot2
#'
#' Internal helper used by .heat_kinship().
#'
#' @param hc An object of class hclust.
#' @param leaf_order Optional leaf labels in the heatmap order.
#'
#' @return A data.frame with x, y, xend, yend.
#'
#' @importFrom stats as.dendrogram
#' @keywords internal
.dendrogram_segments <- function(hc, leaf_order = NULL) {
  dend <- as.dendrogram(hc)

  if (is.null(leaf_order)) {
    leaf_order <- labels(dend)
  }

  leaf_x <- setNames(seq_along(leaf_order), leaf_order)

  segments <- list()

  walk <- function(node) {
    if (is.leaf(node)) {
      lab <- attr(node, "label")
      return(list(x = leaf_x[[lab]], y = 0))
    }

    children <- lapply(node, walk)
    node_height <- attr(node, "height")

    child_x <- vapply(children, `[[`, numeric(1), "x")
    child_y <- vapply(children, `[[`, numeric(1), "y")

    segments[[length(segments) + 1L]] <<- data.frame(
      x = min(child_x),
      y = node_height,
      xend = max(child_x),
      yend = node_height
    )

    for (i in seq_along(children)) {
      segments[[length(segments) + 1L]] <<- data.frame(
        x = child_x[i],
        y = child_y[i],
        xend = child_x[i],
        yend = node_height
      )
    }

    list(x = mean(range(child_x)), y = node_height)
  }

  walk(dend)

  out <- do.call(rbind, segments)
  rownames(out) <- NULL
  out
}


#' Automatic label size based on matrix dimension
#'
#' Internal helper used by .heat_kinship().
#'
#' @param n Number of rows/columns.
#' @param requested_size User-requested label size.
#'
#' @return Numeric ggplot2 text size.
#'
#' @keywords internal
.auto_label_size <- function(n, requested_size = 2.5) {
  multiplier <- requested_size / 2.5

  size <- if (n <= 10) {
    11
  } else if (n <= 15) {
    9.5
  } else if (n <= 20) {
    8
  } else if (n <= 30) {
    6.5
  } else if (n <= 40) {
    5.5
  } else if (n <= 50) {
    4.5
  } else {
    3.5
  }

  size * multiplier
}


#' Kinship heatmaps (adapted from superheat)
#'
#' Implements a subset of superheat behavior used by kinship.heatmap().
#'
#' @param K A numeric matrix.
#' @param dendrogram Logical. Draw a column dendrogram.
#' @param dist.method Distance method passed to dist().
#' @param clustering.method Clustering method. Supports "hierarchical" and "kmeans".
#' @param row.label Logical. Draw row labels.
#' @param col.label Logical. Draw column labels.
#'
#' @return Invisibly returns a list, with plot,
#' membership.cols, membership.rows, order.rows, order.cols, and heat.pal.values.
#'
#' @details Credits to the superheat package by Barter and Yu (2020).
#'
#' @keywords internal

.heat_kinship <- function(
  K,
  dendrogram = TRUE,
  clustering.method = c("hierarchical", "kmeans"),
  dist.method = c(
    "euclidean",
    "maximum",
    "manhattan",
    "canberra",
    "binary",
    "minkowski"
  ),
  row.label = TRUE,
  col.label = FALSE
) {
  clustering.method <- match.arg(clustering.method)
  dist.method <- match.arg(dist.method)

  if (clustering.method == "kmeans")
    warning("kmeans clustering is not implemented, defaulting to hierarchical")

  X <- K

  left.label <- if (isTRUE(row.label)) {
    "variable"
  } else {
    "none"
  }

  bottom.label <- if (isTRUE(col.label)) {
    "variable"
  } else {
    "none"
  }

  bottom.label.text.angle <- 90
  bottom.label.text.size <- 2.5
  legend.text.size <- 10
  legend.width <- 2.0
  left.label.text.size <- 2.5

  if (!is.matrix(X)) {
    stop("K should be a matrix.")
  }

  if (!is.numeric(X)) {
    stop("K should be a numeric matrix.")
  }

  if (any(is.na(X))) {
    stop("K contains missing values.")
  }

  n_rows <- nrow(X)
  n_cols <- ncol(X)

  if (n_rows < 2L || n_cols < 2L) {
    stop("X must contain at least two rows and two columns.")
  }

  row_names <- rownames(X)
  col_names <- colnames(X)
  same_axis_order <- n_rows == n_cols && identical(row_names, col_names)

  hc_rows <- NULL
  hc_cols <- NULL

  order.rows <- seq_len(n_rows)
  order.cols <- seq_len(n_cols)

  membership.rows <- seq_len(n_rows)
  membership.cols <- seq_len(n_cols)

  # The old kinship.heatmap() wrapper always asked superheat for pretty
  # row/column ordering, but did not supply cluster counts or memberships.
  # In that call path, clustering.method is accepted but hclust() still
  # determines the heatmap order.
  if (isTRUE(same_axis_order)) {
    hc_rows <- hclust(dist(X, method = dist.method))
    hc_cols <- hc_rows
    order.rows <- hc_rows$order
    order.cols <- order.rows

    # If axis order is not the same, create one ordering vec for each axis.
  } else {
    hc_rows <- hclust(dist(X, method = dist.method))
    order.rows <- hc_rows$order

    hc_cols <- hclust(dist(t(X), method = dist.method))
    order.cols <- hc_cols$order
  }

  # Reorder based on the clustering.
  X_ord <- X[order.rows, order.cols, drop = FALSE]

  # Adjust row/col names (to the same order as X_ord).
  if (!is.null(row_names)) {
    rownames(X_ord) <- row_names[order.rows]
  }

  if (!is.null(col_names)) {
    colnames(X_ord) <- col_names[order.cols]
  }

  row_labels <- rownames(X_ord)
  col_labels <- colnames(X_ord)

  # Fallback labels.
  if (is.null(row_labels)) {
    row_labels <- as.character(seq_len(n_rows))
    rownames(X_ord) <- row_labels
  }

  # Fallback labels.
  if (is.null(col_labels)) {
    col_labels <- as.character(seq_len(n_cols))
    colnames(X_ord) <- col_labels
  }

  # Move from matrix to long format for plotting.
  heat_df <- as.data.frame(as.table(X_ord), stringsAsFactors = FALSE)
  names(heat_df) <- c("row", "col", "value")

  # Create dimension indices.
  heat_df$x <- match(heat_df$col, col_labels)
  heat_df$row_index <- match(heat_df$row, row_labels)

  # Match superheat orientation: the first ordered row is at the bottom,
  # so the matrix diagonal runs from bottom-left to top-right.
  heat_df$y <- heat_df$row_index

  label_n <- max(n_rows, n_cols)
  row_label_size <- .auto_label_size(label_n, left.label.text.size)
  col_label_size <- .auto_label_size(label_n, bottom.label.text.size)

  # Get pallet from quantiles.
  heat_pal <- c(
    "#440154FF",
    "#3B528BFF",
    "#21908CFF",
    "#5DC863FF",
    "#FDE725FF"
  )
  heat_pal_quantiles <- quantile(
    heat_df$value,
    probs = seq(0, 1, length.out = length(heat_pal)),
    na.rm = TRUE
  )
  heat_pal_range <- range(heat_pal_quantiles, finite = TRUE)
  heat_pal_values <- if (
    length(heat_pal_range) == 2L &&
      is.finite(diff(heat_pal_range)) &&
      diff(heat_pal_range) > 0
  ) {
    (heat_pal_quantiles - heat_pal_range[1]) / diff(heat_pal_range)
  } else {
    seq(0, 1, length.out = length(heat_pal))
  }

  # Get breaks for legend.
  fill_limits <- range(X, finite = TRUE)
  legend_breaks <- if (is.finite(diff(fill_limits)) && diff(fill_limits) > 0) {
    signif(
      as.vector(
        quantile(
          seq(fill_limits[1], fill_limits[2], length.out = 100),
          na.rm = TRUE
        )
      ),
      1
    )
  } else {
    fill_limits[1]
  }

  # Bar styling.
  legend_bar_width <- unit(min(0.95, 0.35 * legend.width), "npc")
  legend_bar_height <- unit(0.25, "cm")

  # Initiate ggplot and add tiles.
  pp <- ggplot() +
    geom_tile(
      data = heat_df,
      aes(x = x, y = y, fill = value),
      colour = "white",
      linewidth = 0.4,
      width = 1,
      height = 1
    ) +
    scale_fill_gradientn(
      colours = heat_pal,
      values = heat_pal_values,
      limits = fill_limits,
      breaks = legend_breaks,
      name = NULL,
      guide = guide_colourbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = legend_bar_width,
        barheight = legend_bar_height,
        ticks = TRUE,
        label.position = "bottom"
      )
    ) +
    scale_x_continuous(
      breaks = seq_len(n_cols),
      labels = if (identical(bottom.label, "variable")) {
        col_labels
      } else {
        rep("", n_cols)
      },
      limits = c(0.5, n_cols + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_len(n_rows),
      labels = if (identical(left.label, "variable")) {
        row_labels
      } else {
        rep("", n_rows)
      },
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 12) +
    # Styling.
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.margin = margin(t = 12, r = 0, b = 0, l = 0),
      legend.text = element_text(size = legend.text.size),
      axis.text.y = element_text(
        size = row_label_size,
        colour = "black",
        margin = margin(r = 8)
      ),
      axis.text.x = element_text(
        size = col_label_size,
        colour = "black",
        angle = bottom.label.text.angle,
        hjust = 1,
        vjust = 0.5,
        margin = margin(t = 8)
      ),
      plot.margin = margin(t = 5, r = 20, b = 5, l = 20)
    )

  # Top buffer.
  y_upper <- n_rows + 0.5

  # Add dendrogram on top is requested.
  if (
    isTRUE(dendrogram) &&
      !is.null(hc_cols)
  ) {
    dend_df <- .dendrogram_segments(hc_cols, leaf_order = col_labels)

    dend_height_tiles <- 1.2 * max(1.1, n_rows * 0.12)
    dend_gap_tiles <- 0.15

    dend_heights <- c(dend_df$y, dend_df$yend)
    dend_heights <- dend_heights[is.finite(dend_heights)]
    max_dend_height <- if (length(dend_heights) > 0L) {
      max(dend_heights)
    } else {
      NA_real_
    }

    if (is.finite(max_dend_height) && max_dend_height > 0) {
      dend_df$y <- n_rows +
        0.5 +
        dend_gap_tiles +
        (dend_df$y / max_dend_height) * dend_height_tiles
      dend_df$yend <- n_rows +
        0.5 +
        dend_gap_tiles +
        (dend_df$yend / max_dend_height) * dend_height_tiles
    } else {
      dend_df$y <- n_rows + 0.5 + dend_gap_tiles
      dend_df$yend <- n_rows + 0.5 + dend_gap_tiles
    }

    # Top buffer adjustment.
    y_upper <- n_rows + 0.5 + dend_gap_tiles + dend_height_tiles

    pp <- pp +
      geom_segment(
        data = dend_df,
        aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend
        ),
        inherit.aes = FALSE,
        colour = "black",
        linewidth = 0.35
      )
  }

  pp <- pp +
    coord_fixed(
      xlim = c(0.5, n_cols + 0.5),
      ylim = c(0.5, y_upper),
      expand = FALSE,
      clip = "off"
    )

  # Remove labels if set to none.
  if (identical(left.label, "none")) {
    pp <- pp +
      theme(axis.text.y = element_blank())
  }

  # Remove labels if set to none.
  if (identical(bottom.label, "none")) {
    pp <- pp +
      theme(axis.text.x = element_blank())
  }

  out <- list(
    plot = pp,
    membership.cols = membership.cols[order.cols],
    membership.rows = membership.rows[order.rows],
    order.rows = order.rows,
    order.cols = order.cols,
    heat.pal.values = heat_pal_values
  )

  print(pp)

  invisible(out)
}


#' Enhanced heatmap plot for a kinship matrix K
#'
#' Generates a heatmap with dendrogram based on a provided kinship matrix.
#' This matrix can be a pedigree relationship matrix \eqn{\boldsymbol{A}}, a
#' genomic relationship matrix \eqn{\boldsymbol{G}} or a hybrid relationship
#' matrix \eqn{\boldsymbol{H}}.
#' Individual names should be assigned to \code{rownames} and \code{colnames}.
#' It sorts individuals according to dendrogram in both columns and rows.
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
#' @examples
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#' G[1:5, 1:5]
#'
#' # Plot a subset of the individuals.
#' kinship.heatmap(K = G[1:10, 1:10], dendrogram = TRUE, row.label = TRUE, col.label = TRUE)
#'
#' @export
kinship.heatmap <- function(
  K = NULL,
  dendrogram = TRUE,
  clustering.method = c("hierarchical", "kmeans"),
  dist.method = c(
    "euclidean",
    "maximum",
    "manhattan",
    "canberra",
    "binary",
    "minkowski"
  ),
  row.label = TRUE,
  col.label = FALSE
) {
  # Check if the class of K is matrix
  if (is.null(K) || !inherits(K, "matrix")) {
    stop("K should be a valid object of class matrix.")
  }

  # Check the rownames/colnames
  if (is.null(rownames(K))) {
    stop("Individual names not assigned to rows of matrix K.")
  }

  if (is.null(colnames(K))) {
    stop("Individual names not assigned to columns of matrix K.")
  }

  if ((identical(rownames(K), colnames(K))) == FALSE) {
    stop("Rownames and colnames of matrix K do not match.")
  }

  # Check if there are missing values
  if (any(is.na(K))) {
    stop("Matrix K contains some missing data.")
  }

  clustering.method <- match.arg(clustering.method)
  dist.method <- match.arg(dist.method)

  pp <- .heat_kinship(
    K = K,
    dendrogram = dendrogram,
    clustering.method = clustering.method,
    dist.method = dist.method,
    row.label = row.label,
    col.label = col.label
  )

  return(pp)
}
