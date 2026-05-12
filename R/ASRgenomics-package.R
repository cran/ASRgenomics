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
#' @importFrom stats prcomp var cov2cor aggregate cor cov na.omit setNames quantile is.leaf
#' @importFrom Matrix nearPD
#' @importFrom AGHmatrix Gmatrix
#' @importFrom crayon blue
#' @importFrom cowplot plot_grid
#' @importFrom ellipse ellipse
#' @importFrom stats dist hclust kmeans prcomp
#' @importFrom ggplot2 ggplot geom_tile geom_segment aes scale_fill_gradientn
#' @importFrom ggplot2 guide_colourbar scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 labs theme_void theme element_text margin coord_fixed
#' @importFrom ggplot2 expansion
#' @importFrom grid unit
#' @importFrom grDevices hcl.colors
## usethis namespace: end
NULL

# Add global variables to avoid NOTES.
utils::globalVariables(
  c(
    "value",
    "Value",
    "density",
    "PC1",
    "PC2",
    "Ind",
    "SNP",
    "MAF",
    "het",
    "Row",
    "Col",
    "n.states",
    "..map",
    "state2",
    "state1",
    "code0",
    "code1A",
    "code1B",
    "code2",
    "row",
    "col",
    "value",
    "x",
    "y",
    "xend",
    "yend"
  )
)
