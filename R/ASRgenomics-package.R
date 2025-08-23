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
  c(
    "value", "Value", "density",
    "PC1", "PC2", "Ind", "SNP", "MAF",
    "het", "Row", "Col", "n.states",
    "..map", "state2", "state1", "code0",
    "code1A", "code1B", "code2"
  )
)

###### Diagnostic of the K matrix message order/content issues

# No individuals filtered by duplicate.thr = 0.95, as none were found.
# There are 49 extreme diagonal values, outside < 0.5 and > 3.2
# There are 0 records of possible duplicates, based on: k(i,j)/sqrt[k(i,i)*k(j,j)] >  0.95
# Cleaning K matrix based on diagonal values and/or duplicates.

# Should be:

# There are 49 extreme diagonal values, outside < 0.5 and > 3.2
# There are 0 records of possible duplicates, based on: k(i,j)/sqrt[k(i,i)*k(j,j)] >  0.95
# Cleaning K matrix based on diagonal values and/or duplicates.
# No individuals filtered by duplicate.thr = 0.95, as none were found.  [REMOVE?]
# A total of 49 individuals were removed based on diagonal thresholds   [Or ADD?]

# And this:

# No individuals filtered out by the diagonal thresholds, as none were found.
# No individuals filtered by duplicate.thr = 0.95, as none were found.
# There are 0 extreme diagonal values, outside < 0.5 and > 3.2
# There are 0 records of possible duplicates, based on: k(i,j)/sqrt[k(i,i)*k(j,j)] >  0.95

# Should be:

# There are 0 extreme diagonal values, outside < 0.5 and > 3.2
# There are 0 records of possible duplicates, based on: k(i,j)/sqrt[k(i,i)*k(j,j)] >  0.95
# No individuals filtered out by the diagonal thresholds, as none were found.   [REMOVE?]
# No individuals filtered by duplicate.thr = 0.95, as none were found.          [REMOVE?]
