#' Check the genomic relationship matrix G against
#' the pedigree relationship matrix A or vice versa
#'
#' Assesses a given genomic relationship matrix \eqn{\boldsymbol{G}} against the
#' pedigree relationship matrix \eqn{\boldsymbol{A}}, or vice versa,
#' to determine the matched and mismatched individuals.
#' If requested, it provides the cleaned versions containing only the matched individuals
#' between both matrices. The user should provide the matrices \eqn{\boldsymbol{G}}and
#' \eqn{\boldsymbol{A}} in full form (\eqn{ng \times ng} and \eqn{na \times na}, respectively).
#' Individual names should be assigned to \code{rownames} and \code{colnames} for both matrices.
#'
#' @param G Input of the genomic relationship matrix \eqn{\boldsymbol{G}} in full form (\eqn{ng \times ng}) (default = \code{NULL}).
#' @param A Input of the pedigree relationship matrix \eqn{\boldsymbol{A}} in full form (\eqn{na \times na}) (default = \code{NULL}).
#' @param clean If \code{TRUE} generates new clean \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{A}}
#' matrices in full form containing only matched individuals (default = \code{TRUE}).
#' @param ord If \code{TRUE} it will order by ascending order of individual names
#' both of the clean \eqn{\boldsymbol{A}} and \eqn{\boldsymbol{G}} matrices (default = \code{TRUE}).
#' @param mism If \code{TRUE} generates two data frames with mismatched individual names
#' from the \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{A}} matrices (default = \code{FALSE}).
#' @param RMdiff If \code{TRUE} it generates the matrix (in lower diagonal row-wise sparse form) of matched
#' observations from both the \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{A}} matrices.
#' This matrix can be used to identify inconsistent values between matched matrices, but it can be very large
#' (default = \code{FALSE}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{Gclean}: a matrix with the portion of \eqn{\boldsymbol{G}} containing only matched individuals.}
#' \item{\code{Aclean}: a matrix with the portion of \eqn{\boldsymbol{A}} containing only matched individuals.}
#' \item{\code{mismG}: a vector containing the names of the individuals from matrix \eqn{\boldsymbol{G}} that are
#' missing in matrix \eqn{\boldsymbol{A}}.}
#' \item{\code{mismA}: a vector containing the names of the individuals from matrix \eqn{\boldsymbol{A}} that are
#' missing in matrix \eqn{\boldsymbol{G}}.}
#' \item{\code{RM}: a data frame with the observations from both the \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{A}}
#' matched matrices, together with their absolute relationship difference.}
#' \item{\code{plotG2A}: scatterplot with the pairing of matched pedigree- against genomic-based
#' relationship values. This graph might take a long to plot with large datasets.}
#' }
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
#' G <- G.matrix(M = M.clean, method = "VanRaden", na.string = "-9")$G
#' G[1:5, 1:5]
#' dim(G)
#'
#' # Match G2A.
#' check <- match.G2A(
#'  A = A, G = G,
#'  clean = TRUE, ord = TRUE, mism = TRUE, RMdiff = TRUE)
#' ls(check)
#' dim(check$Aclean)
#' dim(check$Gclean)
#' check$Aclean[1:5, 1:5]
#' check$Gclean[1:5, 1:5]
#' head(check$mismG)
#' head(check$mismA)
#' check$plotG2A
#' head(check$RM)
#' }
#'

match.G2A <- function(A = NULL, G = NULL, clean = TRUE, ord = TRUE,
                      mism = FALSE, RMdiff = FALSE, message = TRUE){

  if (is.null(A) || !inherits(A, "matrix")) {
    stop("A should be a valid object of class matrix.")
  }
  if (is.null(G) || !inherits(G, "matrix")) {
    stop("G should be a valid object of class matrix.")
  }
  # Check the rownames/colnames
  if (is.null(rownames(A))){
    stop("Individual names not assigned to rows of matrix A.")
  }
  if (is.null(colnames(A))){
    stop('Individual names not assigned to columns of matrix A.')
  }
  if ((identical(rownames(A), colnames(A))) == FALSE){
    stop("Rownames and colnames of matrix A do not match.")
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

  # Check for consistency between A and G
  Aind <- row.names(A)
  Gind <- row.names(G)

  if (all(Aind %in% Gind) & message){
    message("All ", ncol(A), " individuals from matrix A match those individuals from matrix G.")
  }
  if (all(Gind %in% Aind) & message){
    message("All ", ncol(G), " individuals from matrix G match those individuals from matrix A.")
  }

  notGenotyped <- which((Aind %in% Gind) == FALSE)
  notPedigree <- which((Gind %in% Aind) == FALSE)

  # If G has different individuals than A, remove these individuals from G
  if (length(notPedigree) > 0) {
    if (message){
      message("Matrix G has ", length(notPedigree), " individuals (out of ", ncol(G), ") NOT present on matrix A.")
    }
    if (clean) {
      Gclean <- G[-notPedigree,-notPedigree]
    } else {
      Gclean <- NULL
    }
    if (mism) {
      rpG <- rownames(G[notPedigree,])
    } else {
      rpG <- NULL
    }
  } else {
    Gclean <- G
    rpG <- NULL
  }

  # If A has different ind than G, remove these ind from A
  if (length(notGenotyped) > 0) {
    if (message){
      message("Matrix A has ", length(notGenotyped), " individuals (out of ", ncol(A), ") NOT present on matrix G.")
    }
    if (clean) {
      Aclean <- A[-notGenotyped, -notGenotyped]
    } else {
      Aclean <- NULL
    }
    if (mism) {
      rpP <- rownames(A[notGenotyped,])
    } else {
      rpP <- NULL
    }
  } else {
    Aclean <- A
    rpP <- NULL
  }

  # Check order of A and G
  if (all(row.names(Aclean) == row.names(Gclean)) == FALSE){
    if (!ord) {
      if (message){
        message("Order of individual names from matched matrices A and G DO NOT agree.")
      }
    }
    if (ord){
      Gclean <- Gclean[order(rownames(Gclean), decreasing=FALSE),
                       order(colnames(Gclean), decreasing=FALSE)]
      Aclean <- Aclean[order(rownames(Aclean), decreasing=FALSE),
                       order(colnames(Aclean), decreasing=FALSE)]
    }
  }

  # TODO this section has to be improved since A.sparse can also be NULL!

  # Sparse form matrices (faster for plots useful for RM matrix)
  G.sparse <- full2sparse(K=Gclean, drop.zero=FALSE)
  A.sparse <- full2sparse(K=Aclean, drop.zero=FALSE)

  # Generating plot of Aclean vs Gclean
  LL <- min(A.sparse[,3], G.sparse[,3])
  UL <- max(A.sparse[,3], G.sparse[,3])

  # Improved version of plot (faster). GG
  p <- ggplot(data.frame(AValue = A.sparse[,3], GValue = G.sparse[,3]),
              aes(x = AValue, y = GValue, color = 'black')) +
         geom_scattermost(xy = cbind(A.sparse[,3], G.sparse[,3]), color = "#0072B2", pointsize = 2) +
         geom_abline(linetype = "dashed") +
         xlim(LL, UL) + ylim(LL, UL)+
         labs(x="Pedigree Relationship (A matrix)",
         y = "Genomic Relationship (G matrix)") +
         theme_classic()

    # RM matrix for diagnostics
  if (isTRUE(RMdiff)) {
    if (!isTRUE(clean)) {
      stop("Option clean must be TRUE to produce RM data frame.")
    }
    RM <- data.frame(A.sparse[,1:2], AValue=round(A.sparse[,3],6), GValue=round(G.sparse[,3],6))
    RM$absdiff <- round(abs(RM$AValue - RM$GValue),6)
    RM$IDRow <- rownames(Aclean)[RM$Row]
    RM$IDCol <- rownames(Aclean)[RM$Col]
    RM$Diag <- 0
    RM$Diag[RM$Row == RM$Col] <- 1
    RM <- RM[c(1,2,6,7,3,4,5,8)]
  } else {
    RM <- NULL
  }

  AValue <- GValue <-  NULL

  return(list(Aclean=Aclean, Gclean=Gclean, mismG=rpG, mismA=rpP, RM=RM, plotG2A=p))


}

# A.sparse <- rbind(A.sparse, A.sparse, A.sparse, A.sparse, A.sparse, A.sparse, A.sparse, A.sparse, A.sparse)
# A.sparse <- rbind(A.sparse, A.sparse)
# G.sparse <- rbind(G.sparse, G.sparse, G.sparse, G.sparse, G.sparse, G.sparse, G.sparse, G.sparse, G.sparse)
# G.sparse <- rbind(G.sparse, G.sparse)

