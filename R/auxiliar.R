#' Function that silences everything (e.g., \code{cat()}, \code{print()}, \code{message()}, ...)
#'
#' @param code code to be silenced.
#'
#' @return None.
#'
#' @keywords internal

silent_ <- function(code) {
  sink(tempfile()) ; on.exit(sink()) ; invisible(force(code))
}

#' Creates a dummy map if not provided
#'
#' @param marker.id vector with names of markers to compose dummy map.
#' @param message logical value indicating whether diagnostic messages should be printed on screen (default = \code{TRUE}).
#'
#' @return Data frame with dummy map. A single chromosome/linkage group is created and marker
#' distances are a sequence from one to the number of markers.
#'
#' @keywords internal

dummy.map_ <- function(marker.id =  NULL, message = TRUE) {

  # Report if required.
  if (message) message("Creating dummy map.")

  # Get dummy map.
  map <- data.frame(marker = marker.id,
                    chrom = 1,
                    pos = seq_along(marker.id))
}


#' Assess condition of the inverse of \strong{K}
#'
#' @param Kinv An inverse relationship matrix in full or sparse form.
#'
#' @return An object of class character with the condition of the matrix:
#' well-conditioned or ill-conditioned.
#'
#' @keywords internal

# TODO check with can be further improved.
Kinv.condition <- function(Kinv){

  if(nrow(Kinv) != ncol(Kinv)){
    Kinv <- sparse2full(Kinv)
  }

  n <- ncol(Kinv)
  CN.1 <- Kinv[1, 2]/sqrt(Kinv[1, 1] * Kinv[1, 1])
  CN.N <- Kinv[(n - 1), n]/sqrt(Kinv[(n - 1), (n - 1)] * Kinv[n, n])
  max_diag <- abs(max(diag(Kinv)))
  max_off_diag <- max(abs(Kinv - diag(Kinv)))
  if (abs(CN.1) > 0.99 | abs(CN.N) > 0.99 | max_diag > 1000 | max_off_diag > 1000) {
    status <- "ill-conditioned"
  }
  else {
    status <- "well-conditioned"
  }

  return(status)
}
