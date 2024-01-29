#' @keywords internal
"_PACKAGE"

#' @aliases ASRgenomics
#'
## usethis namespace: start
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
## usethis namespace: end
NULL

# Add global variables to avoid NOTES.
utils::globalVariables(
  c("value", "Value", "density",
    "PC1", "PC2", "Ind", "SNP", "MAF",
    "het", "Row", "Col", "n.states",
    "..map", "state2", "state1", "code0",
    "code1A", "code1B", "code2"))

