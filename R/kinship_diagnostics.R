#' Reports summary statistics, plots and filter options for a given kinship matrix K
#'
#' It reports summary statistics, plots and allows for some filter options
#' for diagonal and off-diagonal elements for a given kinship matrix.
#' The input matrix can be a pedigree-based
#' relationship matrix \eqn{\boldsymbol{A}}, a genomic relationship matrix \eqn{\boldsymbol{G}} or a
#' hybrid relationship matrix \eqn{\boldsymbol{H}}.
#' Individual names should be assigned to \code{rownames} and \code{colnames}.
#'
#' @param K Input of a kinship matrix in full format (\eqn{n \times n}) (default = \code{NULL}).
#' @param diagonal.thr.large A threshold value to flag large diagonal values (default = \code{1.2}).
#' @param diagonal.thr.small A threshold value to flag small diagonal values (default = \code{0.8}).
#' @param duplicate.thr A threshold value to flag possible duplicates. Any calculation larger than the
#' threshold based on
#' \eqn{\boldsymbol{k}_{i,i}\mathbin{/}\sqrt{\boldsymbol{k}_{i,i} \times \boldsymbol{k}_{j,j}}}
#' is identified as a duplicate (default = \code{0.95}).
#' @param clean.diagonal If \code{TRUE} returns a kinship matrix filtered by values smaller than
#' \code{diagonal.thr.large} and larger than \code{diagonal.thr.small} (default = \code{FALSE}).
#' @param clean.duplicate If \code{TRUE} return a kinship matrix without the flagged duplicate individuals.
#' All individuals involved are removed (default = \code{FALSE}).
#' @param plots If \code{TRUE} generates graphical output of the diagonal and off-diagonal
#' values of the kinship matrix (default = \code{TRUE}).
#' @param sample.plot A numeric value between 0 and 1 indicating the proportion
#' of the data points to be sampled for fast plotting of off-diagonal values.
#' Note that for proportions other than 1, the method is not exact and low
#' proportions are only recommended for large kinship matrices (default = \code{1}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{list.diagonal}: a data frame with the list of flagged large or small diagonal values.}
#' \item{\code{list.duplicate}: a data frame with the list of possible duplicates.}
#' \item{\code{clean.kinship}: output of kinship matrix filtered without the flagged diagonal
#'  and/or duplicate individuals.}
#' \item{\code{plot.diagonal}: histogram with the distribution of diagonal values from the kinship matrix.}
#' \item{\code{plot.offdiag}: histogram with the distribution of off-diagonal values from kinship matrix.}
#' }
#'
#' @export
#'
#' @examples
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#'
#' # Diagnose G.
#' G_summary <- kinship.diagnostics(
#'  K = G,
#'  diagonal.thr.large = 1.3, diagonal.thr.small = 0.7, clean.diagonal = TRUE,
#'  duplicate.thr = 0.8, clean.duplicate = TRUE,
#'  sample.plot = 0.50)
#' ls(G_summary)
#' dim(G_summary$clean.kinship)
#' G_summary$clean.kinship[1:5, 1:5]
#' G_summary$list.duplicate
#' G_summary$list.diagonal
#' G_summary$plot.diag
#' G_summary$plot.offdiag
#'

kinship.diagnostics <- function(K = NULL,
                                diagonal.thr.large = 1.2, diagonal.thr.small = 0.8,
                                duplicate.thr = 0.95, clean.diagonal = FALSE,
                                clean.duplicate = FALSE,
                                plots = TRUE, sample.plot = 1, message = TRUE){

  # Check if the class of K is matrix.
  if (is.null(K) || !inherits(K, "matrix")) {
    stop("K should be a valid object of class matrix")
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
  # Check on other input
  if (duplicate.thr < 0 | duplicate.thr > 1) {
    stop("Specification of duplicate.thr must be between 0 and 1.")
  }
  if (diagonal.thr.large < diagonal.thr.small) {
    stop("Value of diagonal.thr.large has to be equal or larger than diagonal.thr.small.")
  }
  if (diagonal.thr.large < 0 | diagonal.thr.small < 0) {
    stop("Values of diagonal.thr.large and diagonal.thr.small have the be both positive.")
  }
  if (sample.plot <= 0 | sample.plot > 1) {
    stop("Values of sample.plot must be between 0 and 1.")
  }

  # Preparing submatrices
  n <- nrow(K)
  #indNames <- rownames(K)
  diagK <- diag(K)  # Vector of diagonal
  #offdiag <- K[lower.tri(K, diag=FALSE)]
  Kcorr <- cov2cor(K)
  K.sparse <- full2sparse(K)
  corrS <- full2sparse(Kcorr)
  K.sparse <- data.frame(K.sparse, Corr=corrS[,3])
  offK <- K.sparse[K.sparse$Row != K.sparse$Col,]
  rm(K.sparse, Kcorr, corrS)

  # SOME STATISTICS
  # Some general statistics and report in 'rp'
  if (message){
    message("Matrix dimension is: ", n, "x", n)
  }

  rank <- try(qr(K)$rank, silent=TRUE)

  # Commented out
  #if (message){
  #  if(class(rank) == "try-error"){
  #    message("Rank cannot be obtained due to missing values!")
  #  } else {
  #    message("Rank of matrix is: ",rank)
  #  }
  #}

  range.diagonal <- c(min=min(diagK, na.rm=TRUE), max=max(diagK, na.rm=TRUE))
  if (message){
    message("Range diagonal values: ", round(range.diagonal[1], 5),
            " to ",round(range.diagonal[2], 5))
  }
  mean.diag <- mean(diagK, na.rm=TRUE)
  if (message){
    message("Mean diagonal values: ", round(mean.diag, 5))
  }
  range.off.diagonal <- c(min=min(offK$Value, na.rm=TRUE), max=max(offK$Value, na.rm=TRUE))
  if (message){
    message("Range off-diagonal values: ", round(range.off.diagonal[1], 5),
            " to ",round(range.off.diagonal[2], 5))
  }
  mean.off.diag <- mean(offK$Value, na.rm=TRUE)
  if (message){
    message("Mean off-diagonal values: ", round(mean.off.diag, 5))
  }

  ##########################
  # DEALING with diagonal

  df.list.diag <- data.frame(
    value = sort(diagK[diagK > diagonal.thr.large
                       | diagK < diagonal.thr.small], decreasing=TRUE))

  # Generating list of flagged potential duplicates
  df.list.duplicate <- offK[offK$Corr > duplicate.thr,]
  df.list.duplicate <- data.frame(df.list.duplicate, Indiv.A=rownames(K)[df.list.duplicate$Row],
                                  Indiv.B=colnames(K)[df.list.duplicate$Col])
  rownames(df.list.duplicate) <- NULL
  df.list.duplicate <- df.list.duplicate[,c(5,6,3,4)]
  df.list.duplicate <- df.list.duplicate[order(df.list.duplicate$Corr, decreasing=TRUE),]
  rownames(df.list.duplicate) <- NULL

  # Generating new K if requested with clean.diagonal.
  if (clean.diagonal){
    if (nrow(df.list.diag) > 0){
      Kclean <- K[-which(rownames(K) %in% row.names(df.list.diag)),
                  -which(rownames(K) %in% row.names(df.list.diag))]
    } else {
      if (message){
        message("No individuals filtered out by the diagonal thresholds, as none were found.")
      }

      # I added this line because it will be NULL if nothing has changed.
      Kclean <- NULL
    }
  } else {
    Kclean <- NULL
  }

  # Generating new K if requested with clean.duplicate.
  if (clean.duplicate){
    if (nrow(df.list.duplicate) > 0){
      idx.offdiag <- unique(c(df.list.duplicate$Indiv.A, df.list.duplicate$Indiv.B))
      # idx.offdiag <- unique(df.list.duplicate$Indiv.A)
      if (is.null(Kclean)) {
        Kclean <- K[-which(rownames(K) %in% idx.offdiag),
                    -which(rownames(K) %in% idx.offdiag)]
      } else {
        Kclean <- Kclean[-which(rownames(Kclean) %in% idx.offdiag),
                         -which(rownames(Kclean) %in% idx.offdiag)]
      }
    } else {
      if (message){
        message("No individuals filtered by duplicate.thr = ", duplicate.thr,", as none were found.")
      }
    }
  }

  # Removing K
  rm(K)

  # SOME STATISTICS
  # Report extreme cases
  count <- length(diagK[diagK > diagonal.thr.large | diagK < diagonal.thr.small])
  if (message){
    message("There are ", count, " extreme diagonal values, outside < ", diagonal.thr.small,
            " and > ", diagonal.thr.large)
  }
  count <- nrow(df.list.duplicate)
  if (message){
    message("There are ", count, " records of possible duplicates, based on: k(i,j)/sqrt[k(i,i)*k(j,j)] >  ", duplicate.thr)
  }

  # Obtain histogram of diagonal ----------------------------------------------------------------

  if (plots){

    # Get data.frame for plotting.
    diagK <- as.data.frame(diagK)
    names(diagK)  <- "value"

    # Generate main plot.
    p1 <- ggplot(diagK, aes(x=value)) +
      geom_histogram(aes(y = after_stat(density)), fill='#0072B2', bins=40) +
      geom_density(alpha=0.3, fill="grey", position='identity') +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            plot.title = element_text(face = "bold")) +
      ggtitle("Diagonal Values")

    # Generate boxplot.
    p2 <- ggplot(aes(value), data = diagK) +
      geom_boxplot() +
      theme_void()

    # Combine plots.
    plot.diag <- plot_grid(p1, p2, rel_heights = c(1, 0.2), ncol = 1, nrow = 2, align = "hv")

    # Obtain histogram of off-diagonal ------------------------------------------------------------

    # Sampling for plot if requested.
    if (sample.plot != 1){
      offK <- offK[sample(
        x = 1:nrow(offK),
        size = floor(sample.plot * nrow(offK)),
        replace = FALSE), , drop = FALSE]
    }

    p1 <- ggplot(offK, aes(x = Value)) +
      geom_histogram(aes(y = after_stat(density)), fill = '#0072B2', bins = 40) +
      geom_density(alpha = 0.3, fill = "grey", position = 'identity') +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            plot.title = element_text(face = "bold")) +
      ggtitle("Off-diagonal Values")

    # Boxplot
    p2 <- ggplot(aes(Value), data = offK) +
      geom_boxplot() +
      theme_void()

    # Combine plots.
    plot.offdiag <- plot_grid(p1, p2, rel_heights = c(1, 0.2), ncol = 1, nrow = 2, align = "hv")

  } else {
    # Nulify plots if not requested.
    plot.diag <- plot.offdiag <- NULL
  }


  # Finalize ------------------------------------------------------------------------------------

  if (nrow(df.list.duplicate) == 0) {df.list.duplicate <- NULL}
  if (nrow(df.list.diag) == 0) {df.list.diag <- NULL}

  return(list(list.diagonal=df.list.diag, list.duplicate=df.list.duplicate,
              clean.kinship=Kclean, plot.diag=plot.diag, plot.offdiag=plot.offdiag))

}
