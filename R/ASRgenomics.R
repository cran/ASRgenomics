#'
#' ASRgenomics: A package with complementary genomic functions
#'
#' @description
#'
#' This package presents a series of molecular and genetic routines in the R
#' environment with the aim of assisting in analytical pipelines before and after
#' use of ASReml-R or another library to perform analyses such as Genomic Selection
#' (GS) or Genome-Wide Association Analyses (GWAS).
#'
#' The main tasks considered are:
#' \itemize{
#' \item Preparing and exploring pedigree, phenotypic and genomic data.
#' \item Calculating and evaluating genomic matrices and their inverse.
#' \item Complementing and expanding results from genomic analyses.
#' }
#'
#' The functions implemented consider aspects such as: filtering SNP data for
#' quality control; assessing a kinship matrix \eqn{\boldsymbol{K}} by
#' reporting diagnostics (statistics and plots); performing Principal
#' Component Analyses (PCA) based on
#' kinship or SNP matrices for understanding population structure;
#' calculating the genomic matrix \eqn{\boldsymbol{G}} and its inverse
#' and assessing their quality; matching pedigree- against genomic-based matrices;
#' tuning up a genomic matrix (bend, blend or align); and obtaining the hybrid
#' matrix \eqn{\boldsymbol{H}} as required with single-step GBLUP (ssGBLUP).
#'
#' For the full manual please visit: \url{https://asreml.kb.vsni.co.uk/wp-content/uploads/sites/3/ASRgenomics_Manual.pdf}
#'
#' @docType package
#' @name ASRgenomics
#'
#' @import data.table
#' @import factoextra
#' @import ggplot2
#' @import scattermore
#' @importFrom utils head
#' @importFrom methods new getFunction
#' @importFrom stats prcomp var cov2cor aggregate cor cov na.omit
#' @importFrom Matrix nearPD
#' @importFrom AGHmatrix Gmatrix
#' @importFrom crayon blue
#' @importFrom cowplot plot_grid
#' @importFrom ellipse ellipse
#' @importFrom superheat superheat

NULL

# \if{html}{
#   \out{<div style="cursor: pointer; text-align: center;" onclick="window.location='https://vsni.co.uk/';">}\figure{VSNi_Logo.png}{options: style = "max-width:10\%;"}\out{</div>}
# }

# Add global variables to avoid NOTES.
utils::globalVariables(
  c("value", "Value", "density",
    "PC1", "PC2", "Ind", "SNP", "MAF",
    "het", "Row", "Col", "n.states",
    "..map", "state2", "state1", "code0",
    "code1A", "code1B", "code2"))
