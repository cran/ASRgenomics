#' Performs a Principal Component Analysis (PCA) based on a kinship matrix K
#'
#' Generates a PCA and summary statistics from a given kinship matrix for
#' population structure. This matrix
#' can be a pedigree-based relationship matrix \eqn{\boldsymbol{A}}, a genomic
#' relationship matrix \eqn{\boldsymbol{G}} or a hybrid relationship matrix
#' \eqn{\boldsymbol{H}}. Individual names should be assigned to \code{rownames} and
#' \code{colnames}. There is additional output such as plots and other data frames
#' to be used on other downstream analyses (such as GWAS).
#'
#' It calls function \code{eigen()} to obtain eigenvalues and later generate the PCA and the
#' \code{factoextra} R package to extract and visualize results.
#'
#' @param K Input of a kinship matrix in full form (\eqn{n \times n}) (default = \code{NULL}).
#' @param scale If \code{TRUE} the PCA analysis will scale the kinship matrix, otherwise
#' it is used in its original scale (default = \code{TRUE}).
#' @param label If \code{TRUE} then includes in output individuals names (default = \code{FALSE}).
#' @param ncp The number of PC dimensions to be shown in the screeplot, and to provide
#' in the output data frame (default = \code{10}).
#' @param groups Specifies a vector of class factor that will be used to define different
#' colors for individuals in the PCA plot. It must be presented in the same order as the individuals
#' in the kinship matrix (default = \code{NULL}).
#' @param ellipses If \code{TRUE}, ellipses will will be drawn around each of the define levels in
#' \code{groups} (default = \code{FALSE}).
#'
#' @return A list with the following four elements:
#' \itemize{
#' \item \code{eigenvalues}: a data frame with the eigenvalues and its variances associated with each dimension
#' including only the first \code{ncp} dimensions.
#' \item \code{pca.scores}: a data frame with scores (rotated observations on the new components) including
#' only the first \code{ncp} dimensions.
#' \item \code{plot.pca}: a scatterplot with the first two-dimensions (PC1 and PC2) and their scores.
#' \item \code{plot.scree}: a barchart with the percentage of variances explained by the \code{ncp} dimensions.
#' }
#'
#' @export
#'
#' @examples
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#' G[1:5, 1:5]
#'
#' # Perform the PCA.
#' G_pca <- kinship.pca(K = G, ncp = 10)
#' ls(G_pca)
#' G_pca$eigenvalues
#' head(G_pca$pca.scores)
#' G_pca$plot.pca
#' G_pca$plot.scree
#'
#' # PCA plot by family (17 groups).
#' grp <- as.factor(pheno.apple$Family)
#' G_pca_grp <- kinship.pca(K = G, groups = grp, label = FALSE, ellipses = FALSE)
#' G_pca_grp$plot.pca
#'

kinship.pca <- function(K=NULL, scale = TRUE, label = FALSE, ncp = 10,
                        groups = NULL, ellipses = FALSE){

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
  # Check if the number of individuals is greater than 3
  if (nrow(K) < 3 | ncol(K) < 3){
      stop("Matrix K needs at least 3 individuals.")
  }
  # Check if the are missing values
  if (any(is.na(K))){
    stop("K matrix contains some missing data.")
  }
  if (ncp < 0 | ncp > nrow(K)) {
    stop("Value ncp must be positive and smaller than the number of rows in matrix K.")
  }

  # Generating the pca
  if (scale) {
    K <- cov2cor(K)
  }

  Peig <- eigen(K) # PCA-eigenvalues
  loadings <- Peig$vectors  # Loadings
  PCAcoord <- as.matrix(K) %*% as.matrix(loadings)  # New PCA coordinates
  colnames(PCAcoord) <- paste('PC',c(1:ncol(K)),sep='')
  colnames(loadings) <- paste('PC',c(1:ncol(K)),sep='')
  rownames(loadings) <- rownames(K)

  # Preparing the prcomp class object
  sdev <- sqrt(Peig$values)
  rotation <- loadings
  x <- PCAcoord
  pca <- list(sdev=sdev, rotation=rotation, center=scale, scale=scale, x=x)
  class(pca) <- "prcomp"

  # Percentage of variances explained by each principal component
  scree_plot <- fviz_eig(pca, addlabels=TRUE, ncp=ncp,
                                     barfill = "#0072B2",
                                     barcolor = "#0072B2",
                                     ggtheme = theme_classic())
  # Extract the eigenvalues/variances of the principal dimensions
  eig_var <- get_eig(pca)

  # Plot PCA
  if (isTRUE(label)) {
    if(is.null(groups)) {
       pca_plot <- fviz_pca_ind(pca, geom=c("point","text"),
                                            repel=TRUE,
                                            col.ind = "#0072B2",
                                            ggtheme = theme_classic())
    } else {
       pca_plot <- fviz_pca_ind(pca,
                        geom = c("point","text"),
                        repel = TRUE,
                        col.ind = groups, # color by groups
                        mean.point = FALSE,
                        legend.title = "Groups",
                        ggtheme = theme_classic())
    }
  }
  if (isFALSE(label)) {
    if(is.null(groups)) {
       pca_plot <- fviz_pca_ind(pca, geom="point",
                                            col.ind = "#0072B2",
                                            ggtheme = theme_classic())
    } else {
       pca_plot <- fviz_pca_ind(pca,
                        geom = "point",
                        col.ind = groups, # color by groups,
                        mean.point = FALSE,
                        legend.title = "Groups",
                        ggtheme = theme_classic())
    }
  }

  # Process ellipses if requested.
  if (!is.null(groups) & ellipses){

    # TODO this can be more memory efficient.

    group.comps <- cbind.data.frame(pca$x[, 1:2], groups)

    # Get centroids.
    centroids <- aggregate(cbind(PC1, PC2) ~ groups, data = group.comps, FUN = mean)

    # Get ellipses.
    ellipses.data  <- do.call(
      rbind,
      lapply(unique(group.comps$groups), function(t) {

        data.frame(
          groups = as.character(t),
          ellipse(
            cov(group.comps[group.comps$groups == t, 1:2]),
            centre = as.matrix(centroids[t, 2:3]), level = 0.95),
          stringsAsFactors=FALSE)
      }
      )
    )

    # Add ellipses to plot.
    pca_plot <- pca_plot +
      geom_path(data = ellipses.data, linewidth = .5, inherit.aes = F,
                aes(x = PC1, y = PC2, color = groups))
  }

  # Scores (rotated X observations on the new components) for ncp components
  scores <- pca$x[,c(1:ncp)]
  eigenvalues <- eig_var[c(1:ncp),]

  return(list(pca.scores=scores, eigenvalues=eigenvalues, plot.scree=scree_plot, plot.pca=pca_plot))

  }
