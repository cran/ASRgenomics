#' Recodes the molecular matrix M for downstream analyses
#'
#' Reads molecular data in format of bi-allelic nucleotide bases (AA,
#' AG, GG, CC, etc.) and recodes them as 0, 1, 2 and \code{NA} to be used in other
#' downstream analyses.
#'
#' @param M A character matrix with SNP data of full form (\eqn{n \times p}),
#' with \eqn{n} individuals and \eqn{p} markers
#' Individual and marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' Data in matrix is coded as AA, AG, GG, CC, etc (default = \code{NULL}).
#' @param recoding A character indicating the recoding option to be performed.
#' Currently, only the nucleotide bases (AA, AG, ...) to allele count is available (\code{"ATGCto012"})
#' (default = \code{"ATGCto012"}).
#' @param map (Optional) A data frame with the map information with \eqn{p} rows.
#' If \code{NULL} a dummy map is generated considering a single chromosome and sequential
#' positions for markers and includes reference allele and alternative allele (default = \code{NULL}).
#' @param marker A character indicating the name of the column in data frame \code{map} with the identification
#' of markers. This is mandatory if \code{map} is provided (default = \code{NULL}).
#' @param ref A character indicating the name of the column in the map containing the reference allele for
#' recoding. If absent, then conversion will be based on the major allele (most frequent).
#' The marker information of a given individual with two of the specified major alleles
#' in \code{ref} will be coded as 2. This is mandatory if \code{map} is provided (default = \code{NULL}).
#' @param alt A character indicating the name of the column in the map containing the alternative allele for
#' recoding. If absent, then it will be inferred from the data. The marker information of a given individual
#' with two of the specified alleles in \code{alt} will be coded as 0 (default = \code{NULL}).
#' @param na.string A character that is interpreted as missing values (default = \code{"NA"}).
#' @param rename.markers If \code{TRUE} marker names (as provided in \strong{M}) will be expanded
#' to store the reference and alternative alleles. For example, from AX-88234566 to AX-88234566_C_A.
#' In the event of unidentified alleles, 0 will be used (default = \code{TRUE}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with the following two elements:
#' \itemize{
#' \item \code{Mrecode}: the molecular matrix \eqn{\boldsymbol{M}} recoded to
#' 0, 1, 2 and \code{NA}.
#' \item \code{mapr}: the data frame with the map information including
#' reference allele and alternative allele.
#' }
#'
#' @export
#'
#' @examples
#' # Create bi-allelic base data set.
#' Mnb <- matrix(c(
#'   "A-",  NA, "GG",   "CC",   "AT",   "CC",   "AA",   "AA",
#'   "AAA", NA, "GG",   "AC",   "AT",   "CG",   "AA",   "AT",
#'   "AA",  NA, "GG",   "CC",   "AA",   "CG",   "AA",   "AA",
#'   "AA",  NA, "GG",   "AA",   "AA",    NA,    "AA",   "AA",
#'   "AT",  NA, "GG",   "AA",   "TT",   "CC",   "AT",   "TT",
#'   "AA",  NA,   NA,   "CC",    NA,    "GG",   "AA",   "AA",
#'   "AA",  NA,   NA,   "CC",   "TT",   "CC",   "AA",   "AT",
#'   "TT",  NA, "GG",   "AA",   "AA",   "CC",   "AA",   "AA"),
#'   ncol = 8, byrow = TRUE, dimnames = list(paste0("ind", 1:8),
#'                                        paste0("m", 1:8)))
#' Mnb
#'
#' # Recode without map (but map is created).
#' Mr <- snp.recode(M = Mnb, na.string = NA)
#' Mr$Mrecode
#' Mr$map
#'
#' # Create map.
#' mapnb <- data.frame(
#'  marker = paste0("m", 1:8),
#'  reference = c("A", "T", "G", "C", "T", "C", "A", "T"),
#'  alternative = c("T", "G", "T", "A", "A", "G", "T", "A")
#'  )
#'  mapnb
#'
#' # Recode with map without alternative allele.
#' Mr <- snp.recode(M = Mnb, map = mapnb, marker = "marker", ref = "reference",
#'            na.string = NA, rename.markers = TRUE)
#' Mr$Mrecode
#' Mr$map
#'
#' # Notice that the alternative allele is in the map as a regular variable,
#' # but in the names it is inferred from data (which might be 0 (missing)).
#'
#' # Recode with map with alternative allele.
#' Mr <- snp.recode(M = Mnb, map = mapnb, marker = "marker",
#'  ref = "reference", alt = "alternative",
#'  na.string = NA, rename.markers = TRUE)
#' Mr$Mrecode
#' Mr$map # Now the alternative is also on the names.
#'
#' # We can also recode without renaming the markers.
#' Mr <- snp.recode(M = Mnb, map = mapnb, marker = "marker", ref = "reference",
#'            na.string = NA, rename.markers = FALSE)
#' Mr$Mrecode
#' Mr$map # Now the alternative is also on the names.
#'

snp.recode <- function(M = NULL, map = NULL, marker = NULL, ref = NULL, alt = NULL,
                       recoding = c("ATGCto012"),
                       na.string = NA, rename.markers = TRUE,
                       message = TRUE){

  # Traps ---------------------------------------------------------------------------------------

  # Check recoding.
  # recoding is just a placeholder for now.
  recoding <- match.arg(recoding)

  # Check class of M.
  check.data_(data_ = "M", class_ = "matrix")

  # Check if class of values is character.
  check.data.mode_(data_ = "M", mode_ = "character")

  # Create map if not provided (same as in pruning).
  if (!is.null(map)) {

    # Check map class.
    check.data_(data_ = "map", class_ = "data.frame")

    # Check map names. Check mandatory variables in map.
    if(is.null(marker)){stop("The \'marker' option must be specified if \'map' is provided.")}
    if(is.null(ref)){stop("The \'ref' option must be specified if \'map' is provided.")}

    # Check if they are present in the map.
    map.name.hit <- c(marker, ref, alt) %in% names(map)
    if (!all(map.name.hit)){
      stop("Value provided to argument \'", c("marker", "ref", "alt")[!map.name.hit],
           "' does not correspond to a variable in \'map'.")
    }

    # Match map and M.
    if (!identical(as.character(map[[marker]]), colnames(M))){
      stop("map[[marker]] and colnames(M) must be identical.")
    }
  }

  # Check if all values are composed of 2 letters.
  if (any(nchar(M) != 2, na.rm = TRUE)) {
    warning("Marker(s) not compliant with bi-allelic coding: ",
            paste0(colnames(M)[ceiling(which(nchar(M) != 2) / nrow(M))], collapse = ", "),
            ".\n  The respective datapoints have been replaced with NA.")
    M[nchar(M) != 2] <- NA
  }

  special.char <- apply(X = M, MARGIN = 2, FUN = function(col){
    any(grepl(pattern = "[[:punct:]]", x = col))
  })

  # Removing eventual special characters (e.g. A-).
  if (any(special.char)){
    warning("Special characters identified in marker(s): ",
            paste0(names(special.char)[special.char], collapse = ", "),
            ".\n  The respective datapoints have been replaced with NA.")
    M[grepl(pattern = "[[:punct:]]", x = M)] <- NA
  }

  # Body ----------------------------------------------------------------------------------------

  # Replace na.string by NA.
  if (!is.na(na.string)) {
    if (na.string == "NA") { na.string <- NA }
  }
  if (!is.na(na.string)) {
    if (message){
      message('A total of ', sum(M %in% na.string),
              " values were identified as missing with the string ",
              na.string, " and were replaced by NA.")
    }
    M[M %in% na.string] <- NA
  }

  # Function to get sorted states of a marker.
  get.states_ <- function(m = NULL){
    sort( # TODO sorting is not really required now.
      unique(
        unlist(
          strsplit(x = m, split = ""))))
  }

  # Identify the states.
  states <- apply(X = M, MARGIN = 2, FUN = get.states_, simplify = FALSE)

  # Initiate main frame.
  reference.frame <- data.table()

  # Get number of states.
  reference.frame[, n.states := sapply(states, length) ]

  # Get markers names.
  reference.frame[, marker := colnames(M)]

  # Check for more than 2 states.
  if (any(reference.frame$n.states > 2)) {
    stop("Markers with more than two allelic states: ",
         paste0(colnames(M)[which(reference.frame$n.states > 2)], collapse = ", "),".")
  }

  # Add states to the frame.
  reference.frame[, c("state1", "state2"):=
                    list(sapply(states, function(m) m[1]),
                         sapply(states, function(m) m[2]))]

  # Add reference state to reference frame.
  if(is.null(ref)) {

    # Logic: the following code calculates the MAF based on the bases.
    # It pastes all together, and the separate all and counts the composing letters.
    # It this is too slow we might have to change a bit.
    # The ifelse trick is required if a marker only has missing data.
    reference.frame[, ref :=
                      apply(X = M, MARGIN = 2, FUN = function(m){
                        cur.ref <- names(
                          which.max(
                            table(
                              strsplit(
                                paste0(na.omit(m), collapse = ""),
                                "")
                              )
                            )
                          )
                        return(
                          ifelse(test = is.null(cur.ref),
                                 yes = NA,
                                 no = cur.ref))
                      })
    ]

  } else {
    reference.frame[, ref := ..map[[ref]]]

    if (!is.null(alt)){
      reference.frame[, alt := ..map[[alt]]]
    }
  }

  # Identify the alternative state (based on the reference state).
  if (is.null(alt)){
    reference.frame[, alt := ifelse(test = {state2 == ref}, yes = state1, no = state2)]
  }

  # Find wrong coding in reference allele.
  # Find wrong references.
  wrong.code <- reference.frame[state1 != ref & state2 != ref, marker]

  if (length(wrong.code) > 0) {
    stop("The provided reference (\'ref') missmatches the allele codings in: ", wrong.code, ".")
  }

  # Find wrong coding in alternative allele.
  if (!is.null(alt)){

    # Find wrong references.
    wrong.code <- reference.frame[state1 != alt & state2 != alt, marker]
    if (length(wrong.code) > 0) {
      stop("The provided reference (\'alt') missmatches the allele codings in: ", wrong.code, ".")
    }
  } ; rm(wrong.code)

  # TODO try replacing on a new matrix.
  # Create combinations for comparisons.
  # NAs are checked in the alternative state for the first 3, and on the reference for the 4.
  reference.frame[!is.na(alt), code0 := paste0(alt, alt)]
  reference.frame[!is.na(alt), code1A := paste0(ref, alt)]
  reference.frame[!is.na(alt), code1B := paste0(alt, ref)]
  reference.frame[!is.na(ref), code2 := paste0(ref, ref)]

  # Rename markers if requested.
  if (rename.markers){
    # reference.frame[, marker := paste0(marker, "_", ref, "_", alt)]
    reference.frame[, marker :=
                      paste0(marker, "_",
                             replace(x = ref, list = is.na(ref), values = "0"), "_",
                             replace(x = alt, list = is.na(alt), values = "0"))]

  }

  # Replace letters with numbers.
  M <- sapply(1:ncol(M), FUN = function(index){
    m <- M[, index]
    tmp.ref <- reference.frame[index,]
    m[m %in% na.omit(tmp.ref[["code0"]])] <- 0
    m[m %in% na.omit(tmp.ref[["code1A"]])] <- 1
    m[m %in% na.omit(tmp.ref[["code1B"]])] <- 1
    m[m %in% na.omit(tmp.ref[["code2"]])] <- 2
    return(m)
  })

  # Reassign names to M.
  colnames(M) <- reference.frame[["marker"]]

  # Transform to numeric.
  mode(M) <- "numeric"

  # Finalize ------------------------------------------------------------------------------------

  # Report.
  if (message) {
    message("Matrix M was recoded from bi-allelic nucleotide bases to numeric.")
  }

  # Prepare ref to export.
  if (is.null(map)){
    map <- dummy.map_(marker.id = reference.frame[["marker"]], message = FALSE)

    # Add reference and alternative alleles to map.
    map$ref <- reference.frame$ref
    map$alt <- reference.frame$alt

  } else {

    # If map is not NULL and rename is requested. Collect names from reference frame.
    if(rename.markers){
      map[[marker]] <- reference.frame[["marker"]]
    }
  }


  # Return the output list.
  return(list(Mrecode = M, map = map))

}
