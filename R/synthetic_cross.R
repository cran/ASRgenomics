#' Generates a molecular matrix M for hypothetical crosses based on the
#' genomic information of the parents
#'
#' This function generates (or imputes) a molecular matrix for offspring
#' from hypothetical crosses based on the genomic information from the parents.
#' This is a common procedure in species such as maize, where only the parents
#' (inbred lines) are genotyped, and this information is used to generate/impute
#' the genotypic data of each of the hybrid offspring.
#' This function can be also used for bulked DNA analyses, in order to obtain an
#' bulked molecular matrix for full-sib individuals were only parents are genotyped.
#'
#'
#' @param M A matrix with marker data of full form (\eqn{n \times p}), with \eqn{n} individuals
#' (mothers and fathers) and \eqn{p} markers.
#' Individual and marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' Data in matrix is coded as 0, 1, 2 (integer or numeric) (default = \code{NULL}).
#' @param ped A data frame with three columns containing only the pedigree of the hypothetical offspring.
#' (not pedigree of parents)
#' It should include the three columns for individual, mother and father (default = \code{NULL}).
#' @param indiv A character indicating the column in \code{ped} data frame containing the identification
#' of the offspring (default = \code{NULL}).
#' @param mother A character indicating the column in \code{ped} data frame containing the identification
#' of the mother (default = \code{NULL}).
#' @param father A character indicating the column in \code{ped} data frame containing the identification
#' of the father (default = \code{NULL}).
#' @param heterozygote.action Indicates the action to take when heterozygotes are found in a marker.
#' Options are: \code{"useNA"}, \code{"exact"}, \code{"fail"}, and \code{"expected"}.
#' See details for more information (default = \code{"useNA"})
#' @param na.action  Indicates the action to take when missing values are found in a marker.
#' Options are: \code{"useNA"} and \code{"expected"}.
#' See details for more information (default = \code{"useNA"}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @details
#' For double-haploids, almost the totality of the markers (except for genotyping errors)
#' will be homozygotic reads. But in many other cases (including recombinant inbred lines)
#' there will be a proportion of heterozygotic reads. In these case, it is
#' very difficult to infer (impute) the exact genotype of a given offspring individual.
#' For example, if parents are 0 (AA) and 1 (AC) then offsprings will differ given this
#' Mendelian sampling. However, different strategies exist to determine the
#' expected value for that specific cross (if required), which are detailed below
#' using the option \code{heterozygote.action}.
#'
#' \itemize{
#'
#'  \item{If \code{heterozygote.action =  "useNA"},
#' the generated offspring will have, for the heterozygote read, an \code{NA},
#' and no markers are removed.
#' Hence, no attempt will be done to impute/estimate this value.
#' }
#'
#'  \item{If \code{heterozygote.action = "exact"},
#' any marker containing one or more heterozygote reads will be removed.
#' Hence, inconsistent markers are fully removed from the \eqn{\boldsymbol{M}} matrix.
#' }
#'
#'  \item{If \code{heterozygote.action = "fail"},
#' function stops and informs of the presence of heterozygote reads.
#' }
#'
#'  \item{If \code{heterozygote.action = "expected"},
#' then an algorithm is implemented, on the heterozygote read to determine its
#' expected value. For example, if parents are 0 and 1, then the expected value
#' (with equal probability) is 0.5. For a cross between two heterozygotes,
#' the expected value is: \eqn{0(1/4) + 1(1/2) + 2(1/4) = 1}. And for a cross
#' between 1 and 2, the expected value is: \eqn{1(1/2) + 2(1/2) = 1.5}}
#' }
#'
#' Missing value require special treatment, and an imputation strategy is detailed
#' below as indicated using the option \code{na.action}.
#'
#' \itemize{
#'
#' \item{If \code{na.action = "useNA"}, if at least one of the parental reads
#' is missing values for a given marker then it will be assigned as missing for
#' the hypothetical cross. Hence, no attempt will be done to impute/estimate
#' this value.}
#'
#' \item{If \code{na.action = "expected"}, then an algorithm is implemented that
#' will impute the expected read of the cross if the genotype of \strong{one of
#' the parents is missing} (\emph{e.g.}, cross between 0 and NA). Calculations
#' are based on parental allelic frequencies \eqn{p} and \eqn{q} for the given
#' marker. The expressions for expected values are detailed below.}
#'
#' \itemize{
#'
#' \item{If the genotype of the non-missing parent read is 0.
#'
#' \eqn{q^2} (probability that the missing parent is 0) x 0 (expected value of the offspring from a 0 x 0 cross: \eqn{0(1/1)}) +
#'
#' \eqn{2pq} (probability that the missing parent is 1) x 0.5 (expected offspring from a 0 x 1 cross: \eqn{0(1/2) + 1(1/2)}) +
#'
#' \eqn{q^2} (probability that the missing parent is 2) x 1 (expected offspring from a 0 x 2 cross: \eqn{1(1/1)})}
#'
#' \item{If the genotype of the non-missing parent read is 1.
#'
#' \eqn{q^2} (probability that the missing parent is 0) x 0.5 (offspring: \eqn{0(1/2) + 1(1/2)}) +
#'
#' \eqn{2pq} (probability that the missing parent is 1) x 1 (offspring: \eqn{0(1/4) + 1(1/2) + 2(1/4)}) +
#'
#' \eqn{q^2} (probability that the missing parent is 2) x 1.5 (offspring: \eqn{1(1/2) + 2(1/2)})}
#'
#' \item{If the genotype of the non-missing parent read is 2.
#'
#' \eqn{q^2} (probability that the missing parent is 0) x 1 (offspring: \eqn{1(1/1)}) +
#'
#' \eqn{2pq} (probability that the missing parent is 1) x 1.5 (offspring: \eqn{1(1/2) + 2(1/2)}) +
#'
#' \eqn{q^2} (probability that the missing parent is 2) x 2 (offspring: \eqn{2(1/1)})}
#' }
#'
#'
#' Similarly, the calculation of the expected read of a cross when \strong{both parents are missing} is
#' also based on population allelic frequencies for the given marker.
#' The expressions for expected values are detailed below.
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{q^2 \times q^2} (probability that both parents are 0) x 0 (expected value of the offspring from a 0 x 0 cross: 0(1/1)) +
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{2 \times (q^2 \times 2pq)} (probability that the first parent is 0 and the second is 1; this requires
#' the multiplication by 2 because it is also possible that the first parent is 1 and the second is 0)
#' x 0.5 (offspring: \eqn{0(1/2) + 1(1/2)}) +
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{2 \times (q^2 \times p^2)} (this could be 0 x 2 or 2 x 0) x 1 (offspring: \eqn{1(1/1)}) +
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{2pq \times 2pq} (both parents are 1) x 1 (offspring: \eqn{0(1/4) + 1(1/2) + 2(1/4)}) +
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{2 \times (2pq \times q2)} (this could be 1 x 2 or 2 x 1) x 1.5 (offspring: \eqn{1(1/2) + 2(1/2)}) +
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{p^2 \times p^2} (both parents are 2) x 2 (offspring: \eqn{2(1/1)})
#'
#' Note that the use of \code{na.action = "expected"} is recommended when
#' a large number of offspring will conform the hybrid cross (such as
#' with bulked DNA analyses) for family groups with reasonable number of individuals.
#'
#' \strong{Warning}. If \code{"expected"} is used for \code{heterozygote.action} or \code{na.action},
#' direct transformation of the molecular data to other codings (\emph{e.g.},
#' dominance matrix coded as \code{c(0,1,0)}) is not recommended.
#' }
#'
#' @return
#' A molecular matrix \eqn{\boldsymbol{M}} containing the genotypes generated/imputed for the
#' hypothetical cross.
#'
#' @export
#' @md
#'
#' @examples
#' # Create dummy pedigree (using first 10 as parents).
#' ped <- data.frame(
#'  male = rownames(geno.apple)[1:5],
#'  female = rownames(geno.apple)[6:10])
#' ped$offs <- paste(ped$male, ped$female, sep = "_")
#' ped
#'
#' # Select portion of M for parents.
#' Mp <- geno.apple[c(ped$male, ped$female), 1:15]
#'
#' # Get genotype of crosses removing markers with heterozygotes.
#' synthetic.cross(
#'  M = Mp, ped = ped,
#'  indiv = "offs", mother = "female", father = "male",
#'  heterozygote.action = "exact",
#'  na.action = "useNA")
#'
#' # Request the synthetic cross to be NA in the respective samples.
#' synthetic.cross(
#'  M = Mp, ped = ped,
#'  indiv = "offs", mother = "female", father = "male",
#'  heterozygote.action = "useNA",
#'  na.action = "useNA")
#'
#' # Get genotype of crosses and use expected values.
#' synthetic.cross(
#'  M = Mp, ped = ped,
#'  indiv = "offs", mother = "female", father = "male",
#'  heterozygote.action = "expected", na.action = "expected")
#'

synthetic.cross <- function(M = NULL, ped = NULL, indiv = NULL, mother = NULL, father = NULL,
                            heterozygote.action = c("useNA", "exact", "fail", "expected"),
                            na.action = c("useNA", "expected"),
                            # n.cores = 1,
                            message = TRUE){

  # Traps ---------------------------------------------------------------------------------------

  # Check na.action.
  na.action <- match.arg(na.action)

  # Check heterozygote.action.
  heterozygote.action <- match.arg(heterozygote.action)

  # Check M class.
  check.data_(data_ = "M", class_ = "matrix")

  # Check ped class.
  check.data_(data_ = "ped", class_ = "data.frame")

  # Check indiv.
  check.args_(data_ = ped, mandatory_ = TRUE, arg_ = indiv,
              class_ = "character", mutate_ = TRUE,
              class.action_ = "message", message_ = TRUE)

  # Check indiv.
  check.args_(data_ = ped, mandatory_ = TRUE, arg_ = mother,
              class_ = "character", mutate_ = TRUE,
              class.action_ = "message", message_ = TRUE)

  # Check indiv.
  check.args_(data_ = ped, mandatory_ = TRUE, arg_ = father,
              class_ = "character", mutate_ = TRUE,
              class.action_ = "message", message_ = TRUE)

  # # Check ped names.
  # ped.name.hit <- c(indiv, mother, father) %in% names(ped)
  # if (!all(ped.name.hit)){
  #   stop("Value provided to argument \'", c("indiv", "mother", "father")[!ped.name.hit],
  #        "' does not correspond to a variable in \'ped'.")
  # }

  # Check heterozygote action.
  heterozygote.action <- match.arg(heterozygote.action)

  # Check na action.
  na.action <- match.arg(na.action)

  # Check message.
  check.logical_(arg_ = "message")

  # Stop if heterozygous action is fail.
  if (heterozygote.action == "fail"){

    # Identify markers with heterozygotes and remove.
    het.markers <- apply(X = M, MARGIN = 2, FUN = function(m) {any(m == 1, na.rm = TRUE)})
    if (any(het.markers)){
      stop("Stop requested as some of the markers have have heterozygous states (1), e.g., ",
           paste0(names(head(het.markers[het.markers], 5)), collapse = ", "), "...")
    }
  }

  # Check if all parents have genotype in M.
  if(!all(c(ped[[mother]], ped[[father]]) %in% rownames(M))){
    stop("Some parents do not have genotypic information in matrix M.")
  }

  # Data manipulation traps ---------------------------------------------------------------------

  # Get unique combinations.
  unique.comb <- paste(ped[[mother]], ped[[father]], sep = "_")

  # Get reciprocals.
  unique.comb.rec <- paste(ped[[father]], ped[[mother]], sep = "_")
  recs <- unique.comb %in% unique.comb.rec

  # Identify duplicates.
  dups <- duplicated(unique.comb)
  if(any(dups)){

    # Report duplicates.
    message("A total of ", sum(dups), " duplicated rows were found in the supplied pedigree.")
    message("Removing \'ped' lines: ", which(dups), ".")

    # Remove duplicates.
    ped <- ped[!dups, ]

  }

  if(any(recs)){

    # Report reciprocals.
    warning("Reciprocals found in the supplied pedigree. These were not removed.")

    # TODO try to remove reciprocals?
  }

  rm(dups, recs)

  # Prepare data --------------------------------------------------------------------------------

  # Collect initial number of markers.
  initial.n.markers <- ncol(M)

  # Check and remove eventual heterozygotes in data if exact method is chosen.
  if (heterozygote.action == "exact"){

    if (message){
      message(blue("\nExact method selected. Eliminating markers containing one or more heterozygotic read."))
    }

    # Identify markers with heterozygotes and remove.
    het.markers <- apply(X = M, MARGIN = 2, FUN = function(m) {any(m == 1, na.rm = TRUE)})
    M <- M[, !het.markers, drop = FALSE]

    if (message){
      message("Total number of dropped markers: ", initial.n.markers - ncol(M))
      message("Total number of remaining markers: ", ncol(M))
    }

    # Stop if no marker left.

    if (ncol(M) == 0){
      stop("No marker were left after removal of heterozygote reads.")
    }
  }

  # Generate hybrid space -----------------------------------------------------------------------

  # Collect marker names (this is necessary in some cases, e.g., 1 marker left).
  marker.ids <- colnames(M)

  # Get range to loop through.
  range.hybrids <- 1:nrow(ped)

  if (message){
    message(blue("\nGenerating in hypothetical crosses genotypic information."))
  }


  # Call function to get offspring genotype.
  get.off.gen <- function(cur.offspring = NULL){

    # Collect genotype of combination.
    M.mother <- M[ped[cur.offspring, mother],]
    M.father <- M[ped[cur.offspring, father],]
    M.offspring <- rbind(M.mother, M.father)

    if (heterozygote.action == "useNA"){
      M.offspring[M.offspring %in% 1] <- NaN
    }

    # Get genotypes with NA.
    if (na.action == "useNA"){
      return(colMeans(M.offspring))
    }

    # Get genotypes with expected values.
    if (na.action == "expected"){

      # Initiate evo object.
      evo <- NULL

      # Get expected means.
      kid <- colMeans(M.offspring)

      # Identify genotypes missing in the current kid
      missing.markers <- is.na(kid) & !is.nan(kid)

      if (any(missing.markers)){

        # Get the frequency of pseudo-major allele 2.
        freq2 <- colMeans(M, na.rm = T)/2

        # Loop across all missing markers from the current kid.
        evo <- sapply(X = which(missing.markers), function(m) {

          # Get current marker.
          marker <- M.offspring[, m]

          # Get pseudo-p.
          f.pseudo.major <- freq2[m]

          # Get pseudo-q.
          f.pseudo.minor <- (1 - f.pseudo.major)

          # Logic: pseudo-q and p are required so we know which frequency to use on
          # the multiplications below. The MAF is not good here because
          # the major allele might not be represented by a 2 in the molecular matrix as
          # we dont know where this data comes from. Using pseudo-q and p, there is no
          # need to know the reference allele! We have to mach the genotype with
          # the possible offspring.

          # If there is one missing parent.
          if(sum(is.na(marker)) == 1){
            par.gen.no.miss <- sum(marker, na.rm = TRUE)

            # If the genotype of the non-missing parent is 0.
            if (par.gen.no.miss == 0){
              evo <-
                # q2
                f.pseudo.minor^2 * 0 +
                # 2pq
                2 * f.pseudo.minor * f.pseudo.major * 0.5 +
                # q2
                f.pseudo.major^2 * 1
            }

            # If the genotype of the non-missing parent is 1.
            if (par.gen.no.miss == 1){

              evo <-
                # q2
                f.pseudo.minor^2 * 0.5 +
                # 2pq
                2 * f.pseudo.minor * f.pseudo.major * 1 +
                # q2
                f.pseudo.major^2 * 1.5
            }

            # If the genotype of the non-missing parent is 2.
            if (par.gen.no.miss == 2){
              evo <-
                # q2
                f.pseudo.minor^2 * 1 +
                # 2pq
                2 * f.pseudo.minor * f.pseudo.major * 1.5 +
                # q2
                f.pseudo.major^2 * 2
            }
          }

          # If there are two missing parents.
          if(sum(is.na(marker)) == 2){
            f.pseudo.major <- freq2[m]

            # All possible combiations of unknown parents.
            evo <-

              # q2 x q2 (0 x 0)
              f.pseudo.minor^2 * f.pseudo.minor^2   * 0 +

              # 2 * q2 x 2pq (this could be 0 x 1 or 1 x 0)
              2 * (f.pseudo.minor^2 * 2 * f.pseudo.minor * f.pseudo.major)   * 0.5+

              # 2 * q2 x p2 (this could be 0 x 2 or 2 x 0)
              2 * (f.pseudo.minor^2 * f.pseudo.major^2)   * 1 +

              # 2pq x 2pq
              2 * f.pseudo.minor * f.pseudo.major * 2 * f.pseudo.minor * f.pseudo.major   *  1 +

              # 2 * 2pq x q2 (this could be 1 x 2 or 2 x 1)
              2 * (2 * f.pseudo.minor * f.pseudo.major * f.pseudo.major^2)   * 1.5 +

              # p2 x p2
              f.pseudo.major^2 * f.pseudo.major^2   * 2

            # Simplification (not sure if works on all cases - probably yes).
            # evo <- f.pseudo.minor^2 * 0 +
            #   2 * f.pseudo.minor * f.pseudo.major * 1 +
            #   f.pseudo.major^2 * 2

          }

          # Return the expected value of the offspring.
          return(evo)
        })


      }

      # Replace na in kid.
      kid[which(missing.markers)] <- evo

      # Return the imputed kid.
      return(kid)
    }
  }

  # Run function.
  # if (n.cores > 1){
  #   M <- mclapply(X = range.hybrids, mc.cores = n.cores, FUN = get.off.gen)
  # } else {
    M <- lapply(X = range.hybrids, FUN = get.off.gen)
  # }

  # Get hybrid matrix.
  M <- do.call(rbind, M)

  # Replace NaN with NA because of heterozygotes.
  M[is.nan(M)] <- NA

  # Add names of hybrids to new matrix.
  rownames(M) <- ped[[indiv]]

  # Add maker ids.
  colnames(M) <- marker.ids

  # Return.
  return(M)
}

