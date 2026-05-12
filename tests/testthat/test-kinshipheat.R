
# Get data ----------------------------------------------------------------------------------------------

G <- G.matrix(M = geno.apple[1:20, 1:50], method = "VanRaden", na.string="NA")$G

# Check kinship plots -----------------------------------------------------------------------------------

test_that("kinship heatmap works", {
  has_dendrogram_layer <- function(plot) {
    built <- ggplot2::ggplot_build(plot)
    any(vapply(
      built$data,
      function(layer) all(c("xend", "yend") %in% names(layer)),
      logical(1)
    ))
  }

  # Plotting normally.
  expect_no_error(
    heat <-
      kinship.heatmap(
        K = G,
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
    )

  # Plotting normally with less things.
  expect_no_error(
    heat <-
      kinship.heatmap(
        K = G,
        dendrogram = FALSE,
        row.label = FALSE,
        col.label = FALSE)
  )

  expect_warning(
    heat_kmeans <-
      kinship.heatmap(
        K = G,
        dendrogram = FALSE,
        clustering.method = "kmeans",
        row.label = TRUE,
        col.label = TRUE),
    "kmeans clustering is not implemented"
  )
  expect_identical(heat_kmeans$order.rows, heat_kmeans$order.cols)
  expect_identical(
    heat_kmeans$order.cols,
    hclust(dist(G, method = "euclidean"))$order
  )

  # The same accepted-but-order-neutral kmeans argument should still draw
  # the requested column dendrogram.
  expect_warning(
    heat_kmeans_dend <-
      kinship.heatmap(
        K = G,
        dendrogram = TRUE,
        clustering.method = "kmeans",
        row.label = TRUE,
        col.label = TRUE),
    "kmeans clustering is not implemented"
  )
  expect_identical(heat_kmeans_dend$order.rows, heat_kmeans_dend$order.cols)
  expect_identical(
    heat_kmeans_dend$order.cols,
    hclust(dist(G, method = "euclidean"))$order
  )
  expect_true(has_dendrogram_layer(heat_kmeans_dend$plot))
})

test_that("traps work", {

  # Using data.frame.
  expect_error(
    heat <-
      kinship.heatmap(
        K = as.data.frame(G),
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
  )

  # No rownames.
  Gwr <- G
  rownames(Gwr) <- c()
  expect_error(
    heat <-
      kinship.heatmap(
        K = Gwr,
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
  )

  # No colnames.
  Gwr <- G
  colnames(Gwr) <- c()
  expect_error(
    heat <-
      kinship.heatmap(
        K = Gwr,
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
  )

  # Diff row/colnames.
  Gwr <- G
  colnames(Gwr)[1] <- 'nil'
  expect_error(
    heat <-
      kinship.heatmap(
        K = Gwr,
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
  )

  # Missing vals.
  Gwr <- G
  Gwr[1,1] <- NA
  expect_error(
    heat <-
      kinship.heatmap(
        K = Gwr,
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE)
  )


})

test_that("internal kinship heatmap helper covers edge branches", {
  heat_fun <- ASRgenomics:::.heat_kinship
  label_size <- ASRgenomics:::.auto_label_size

  expect_equal(label_size(10), 11)
  expect_equal(label_size(15), 9.5)
  expect_equal(label_size(20), 8)
  expect_equal(label_size(30), 6.5)
  expect_equal(label_size(40), 5.5)
  expect_equal(label_size(50), 4.5)
  expect_equal(label_size(51), 3.5)

  expect_error(
    heat_fun(as.data.frame(G)),
    "K should be a matrix"
  )

  G_char <- matrix(letters[1:4], 2, 2)
  expect_error(
    heat_fun(G_char),
    "K should be a numeric matrix"
  )

  G_na <- G[1:2, 1:2]
  G_na[1, 1] <- NA
  expect_error(
    heat_fun(G_na),
    "K contains missing values"
  )

  expect_error(
    heat_fun(matrix(1, 1, 1)),
    "at least two rows and two columns"
  )

  # Non-square matrices use separate row and column ordering branches.
  G_rect <- matrix(
    c(
      1.00, 0.20, 0.30, 0.40,
      0.15, 1.10, 0.35, 0.45,
      0.25, 0.30, 1.20, 0.55
    ),
    nrow = 3,
    byrow = TRUE
  )

  expect_no_error(
    heat_rect_hier <-
      heat_fun(
        K = G_rect,
        dendrogram = FALSE,
        row.label = TRUE,
        col.label = TRUE
      )
  )
  expect_length(heat_rect_hier$order.rows, nrow(G_rect))
  expect_length(heat_rect_hier$order.cols, ncol(G_rect))
  expect_false(anyNA(ggplot2::ggplot_build(heat_rect_hier$plot)$data[[1]]$x))
  expect_false(anyNA(ggplot2::ggplot_build(heat_rect_hier$plot)$data[[1]]$y))
  heat_pal_quantiles <- quantile(
    as.vector(G_rect),
    probs = seq(0, 1, length.out = 5),
    na.rm = TRUE
  )
  expect_equal(
    as.numeric(heat_rect_hier$heat.pal.values),
    as.numeric(
      (heat_pal_quantiles - min(heat_pal_quantiles)) /
        diff(range(heat_pal_quantiles))
    )
  )

  set.seed(123)
  expect_warning(
    heat_rect_kmeans <-
      heat_fun(
        K = G_rect,
        dendrogram = TRUE,
        clustering.method = "kmeans",
        row.label = FALSE,
        col.label = FALSE
      ),
    "kmeans clustering is not implemented"
  )
  expect_length(heat_rect_kmeans$order.rows, nrow(G_rect))
  expect_length(heat_rect_kmeans$order.cols, ncol(G_rect))
  expect_identical(
    heat_rect_kmeans$order.rows,
    hclust(dist(G_rect, method = "euclidean"))$order
  )
  expect_identical(
    heat_rect_kmeans$order.cols,
    hclust(dist(t(G_rect), method = "euclidean"))$order
  )
  expect_true(any(vapply(
    ggplot2::ggplot_build(heat_rect_kmeans$plot)$data,
    function(layer) all(c("xend", "yend") %in% names(layer)),
    logical(1)
  )))

  expect_warning(
    heat_rect_kmeans_plain <-
      heat_fun(
        K = G_rect,
        dendrogram = FALSE,
        clustering.method = "kmeans",
        row.label = FALSE,
        col.label = FALSE
      ),
    "kmeans clustering is not implemented"
  )
  expect_identical(
    heat_rect_kmeans_plain$order.rows,
    hclust(dist(G_rect, method = "euclidean"))$order
  )
  expect_identical(
    heat_rect_kmeans_plain$order.cols,
    hclust(dist(t(G_rect), method = "euclidean"))$order
  )
  expect_false(any(vapply(
    ggplot2::ggplot_build(heat_rect_kmeans_plain$plot)$data,
    function(layer) all(c("xend", "yend") %in% names(layer)),
    logical(1)
  )))

  # A zero-distance dendrogram covers the zero-height dendrogram scaling branch.
  G_zero <- matrix(1, 3, 3)
  rownames(G_zero) <- colnames(G_zero) <- paste0("id", seq_len(3))
  expect_no_error(
    heat_zero <-
      heat_fun(
        K = G_zero,
        dendrogram = TRUE,
        row.label = TRUE,
        col.label = TRUE
      )
  )
  expect_s3_class(heat_zero$plot, "ggplot")
  expect_equal(heat_zero$heat.pal.values, seq(0, 1, length.out = 5))
})
