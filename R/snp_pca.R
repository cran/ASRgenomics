#' Performs a Principal Component Analysis (PCA) based on a molecular matrix M
#'
#' Generates a PCA and summary statistics from a given molecular matrix
#' for population structure. Matrix
#' provided is of full form (\eqn{n \times p}), with n individuals and p markers. Individual and
#' marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as 0, 1, 2 (integers or decimal numbers). Missing values are
#' not accepted and these need to be imputed (see function \code{qc.filtering()}
#' for implementing mean imputation). There is additional output such as plots and
#' other data frames
#' to be used on other downstream analyses (such as GWAS).
#'
#' It calls function \code{prcomp()} to generate the PCA and the
#' \code{factoextra} R package to extract and visualize results.
#' Methodology uses normalized allele frequencies as proposed by Patterson \emph{et al.} (2006).
#'
#' @param M A matrix with SNP data of full form (\eqn{n \times p}), with \eqn{n}
#' individuals and \eqn{p} markers (default = \code{NULL}).
#' @param label If \code{TRUE} then includes in output individuals names (default = \code{FALSE}).
#' @param ncp The number of PC dimensions to be shown in the screeplot, and to provide
#' in the output data frame (default = \code{10}).
#' @param groups Specifies a vector of class factor that will be used to define different
#' colors for individuals in the PCA plot. It must be presented in the same order as the individuals
#' in the molecular \eqn{\boldsymbol{M}} matrix (default = \code{NULL}).
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
#' @references
#' Patterson N., Price A.L., and Reich, D. 2006. Population structure and eigenanalysis.
#' PLoS Genet 2(12):e190. doi:10.1371/journal.pgen.0020190
#'
#' @export
#'
#' @examples
#' # Perform the PCA.
#' SNP_pca <- snp.pca(M = geno.apple, ncp = 10)
#' ls(SNP_pca)
#' SNP_pca$eigenvalues
#' head(SNP_pca$pca.scores)
#' SNP_pca$plot.pca
#' SNP_pca$plot.scree
#'
#' # PCA plot by family (17 groups).
#' grp <- as.factor(pheno.apple$Family)
#' SNP_pca_grp <- snp.pca(M = geno.apple, groups = grp, label = FALSE)
#' SNP_pca_grp$plot.pca
#'


snp.pca <- function(M = NULL, label = FALSE, ncp = 10,
                    groups = NULL, ellipses = FALSE){

  # Check if the class of M is matrix
  if (is.null(M) || !inherits(M, "matrix")) {
    stop("M should be a valid object of class matrix.")
  }
  if (is.null(colnames(M))){
    stop("Marker names not assigned to columns of matrix M.")
  }
  if (is.null(rownames(M))){
    stop("Individuals names not assigned to rows of matrix M.")
  }
  # Check if the are missing values
  if (any(is.na(M))){
    stop("M matrix contains some missing data, consider performing some imputation.")
  }
  if (ncp < 0 | ncp > nrow(M)) {
    stop("Value ncp must be positive and smaller than the number of rows in matrix M.")
  }

  ## PCA by Petterson et al. (2006)
  M.mean <- colMeans(M)/2
  M.scale <- sqrt(M.mean * (1 - M.mean))
  M.norm <- matrix(NA, nrow = nrow(M), ncol = ncol(M))
  for (i in 1:ncol(M)) {  # Done by SNP column
    M.norm[,i] <- (M[,i]/2 - M.mean[i]) / M.scale[i]
  }

  # Pass name of individuals along.
  rownames(M.norm) <- rownames(M)
  colnames(M.norm) <- colnames(M)

  # Generating the pca
  pca <- prcomp(var(t(M.norm)), scale.=FALSE)  # Original it takes a lont time

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

    # TODO this can be more memory efficient (the data frame is being recreated).

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


  # Scores (rotated X observations on the new components) for ncp components.
  scores <- pca$x[,c(1:ncp)]
  eigenvalues <- eig_var[c(1:ncp),]

  return(list(pca.scores=scores, eigenvalues=eigenvalues, plot.scree=scree_plot, plot.pca=pca_plot))

}
