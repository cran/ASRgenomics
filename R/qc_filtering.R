#' Quality control filtering of molecular matrix M for downstream analyses
#'
#' Reads molecular data in the format 0, 1, 2 and performs some basic quality control
#' filters and simple imputation.
#' Matrix provided is of the full form (\eqn{n \times p}), with \eqn{n} individuals and \eqn{p} markers.
#' Individual and marker names are assigned to \code{rownames} and \code{colnames},
#' respectively. Filtering can be done with the some of the following options by
#' specifying thresholds for:
#' missing values on individuals, missing values on markers, minor allele frequency,
#' inbreeding Fis value (of markers), and observed heterozygosity (of markers).
#' String used for identifying missing values can be specified.
#' If requested, missing values will be imputed based on the mean of each SNP.
#'
#' @param M A matrix with SNP data of full form (\eqn{n \times p}), with \eqn{n} individuals and \eqn{p} markers
#' Individual and marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' Data in matrix is coded as 0, 1, 2 (integer or numeric) (default = \code{NULL}).
#' @param base If \code{TRUE} matrix \eqn{\boldsymbol{M}} is considered as bi-allele SNP data format (character)
#' and the SNPs are recoded to numerical values before performing the quality control filters
#' (default = \code{FALSE}) (currently deprecated).
#' @param na.string A character that will be interpreted as \code{NA} values (default = \code{"NA"}).
#' @param map (Optional) A data frame with the map information with \eqn{p} rows (default = \code{NULL}).
#' @param marker A character indicating the name of the column in data frame \code{map} with the identification
#' of markers. This is mandatory if \code{map} is provided (default = \code{NULL}).
#' @param chrom A character indicating the name of the column in data frame \code{map} with the identification
#' of chromosomes (default = \code{NULL}).
#' @param pos A character indicating the name of the column in data frame \code{map} with the identification
#' of marker positions (default = \code{NULL}).
#' @param ref A character indicating the name of the column in the map containing the reference allele for
#' recoding. If absent, then conversion will be based on the major allele (most frequent).
#' The marker information of a given individuals with two of the specified major alleles
#' in \code{ref} will be coded as 2 (default = \code{NULL}).
#' @param marker.callrate A numerical value between 0 and 1 used to remove SNPs with a rate
#' of missing values equal or larger than this value (default = 1, \emph{i.e.} no removing).
#' @param ind.callrate A numerical value between 0 and 1 used to remove individuals with a
#' rate of missing values equal or larger than this value (default = 1, \emph{i.e.} no removing).
#' @param maf A numerical value between 0 and 1 used to remove SNPs with a Minor Allele Frequency
#' (MAF) below this value (default = 0, \emph{i.e.} no removing).
#' @param heterozygosity A numeric value indicating the maximum value of accepted observed heterozygosity (Ho)
#' (default = 1, \emph{i.e.} no removing).
#' @param Fis A numeric value indicating the maximum value of accepted inbreeding (Fis) following
#' the equation \eqn{|1 - (Ho/He)|} (default = 1, \emph{i.e.} no removing).
#' @param impute If \code{TRUE} imputation of missing values is done using the mean of each SNP
#' (default = \code{FALSE}).
#' @param Mrecode If \code{TRUE} it provides the recoded \eqn{\boldsymbol{M}} matrix from the bi-allelic to numeric SNP
#' (default = \code{FALSE}) (currently deprecated).
#' @param plots If \code{TRUE} generates graphical output of the quality control based on the
#' original input matrix (default = \code{TRUE}).
#' @param digits Set up the number of digits used to round the output matrix (default = 2).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{M.clean}: the cleaned \eqn{\boldsymbol{M}} matrix after the quality control filters have been applied.
#' \item \code{map}: if provided, a cleaned \code{map} data frame after the quality control filters have been applied.
#' \item \code{plot.missing.ind}: a plot of missing data per individual (original marker matrix).
#' \item \code{plot.missing.SNP}: a plot of missing data per SNP (original marker matrix).
#' \item \code{plot.heteroz}: a plot of observed heterozygocity per SNP (original marker matrix).
#' \item \code{plot.Fis}: a plot of Fis per SNP (original marker matrix).
#' \item \code{plot.maf}: a plot of the minor allele frequency (original marker matrix).
#' }
#'
#' @md
#' @details
#' \strong{Warning}: The arguments \code{base}, \code{ref}, and \code{Mrecode}
#' currently are deprecated and will
#' be removed on the next version of \code{ASRgenomics}.
#' Use function \link{snp.recode} to recode the matrix prior to using \code{qc.filtering}.
#'
#' The filtering process is carried out as expressed in the following simplified pseudo-code
#' that consists on a loop repeated twice:
#'
#' \strong{for i in 1 to 2}
#'
#' &nbsp; &nbsp; Filter markers based on call rate.
#'
#' &nbsp; &nbsp; Filter individuals based on call rate.
#'
#' &nbsp; &nbsp; Filter markers based on minor allele frequency.
#'
#' &nbsp; &nbsp; Filter markers based on observed heterozygosity.
#'
#' &nbsp; &nbsp; Filter markers based on inbreeding.
#'
#' \strong{end for}
#'
#' @export
#'
#' @examples
#' # Example: Pine dataset (coded as 0,1,2 with missing as -9).
#'
#' M.clean <- qc.filtering(
#'  M = geno.pine926,
#'  maf = 0.05,
#'  marker.callrate = 0.9, ind.callrate = 0.9,
#'  heterozygosity = 0.9, Fis = 0.6,
#'  na.string = "-9")
#' ls(M.clean)
#' M.clean$M.clean[1:5, 1:5]
#' dim(M.clean$M.clean)
#' head(M.clean$map)
#' M.clean$plot.maf
#' M.clean$plot.missing.ind
#' M.clean$plot.missing.SNP
#' M.clean$plot.heteroz
#' M.clean$plot.Fis
#'
#' \donttest{
#' # Example: Salmon dataset (coded as 0,1,2 with missing as NA).
#'
#' M.clean <- qc.filtering(
#'  M = geno.salmon,
#'  maf = 0.02,
#'  marker.callrate = 0.10, ind.callrate = 0.20,
#'  heterozygosity = 0.9, Fis = 0.4)
#' M.clean$M.clean[1:5, 1:5]
#' dim(M.clean$M.clean)
#' head(M.clean$map)
#' M.clean$plot.maf
#' M.clean$plot.missing.ind
#' M.clean$plot.missing.SNP
#' M.clean$plot.heteroz
#' M.clean$plot.Fis
#' }
#'


qc.filtering <- function(M = NULL, base = FALSE, na.string = NA,
                         map = NULL, marker = NULL, chrom = NULL, pos = NULL, ref = NULL,
                         marker.callrate = 1, ind.callrate = 1, maf = 0,
                         heterozygosity = 1, Fis = 1,
                         impute = FALSE, Mrecode = FALSE,
                         plots = TRUE, digits = 2, message = TRUE) {

  # Deprecation traps ---------------------------------------------------------------------------

  if (!is.null(ref) | base | Mrecode){
    stop("The recoding has been deprecated in \'qc.filtering()', please use \'snp.recode()' to perform this task.")
  }

  # Traps ---------------------------------------------------------------------------------------

  # Check if the class of M is matrix.
  if (is.null(M) || !inherits(M, "matrix"))
    stop("M should be a valid object of class matrix.")

  if (is.null(colnames(M)))
    stop("Marker names not assigned to columns of matrix M.")

  if (is.null(rownames(M)))
    stop("Individuals names not assigned to rows of matrix M.")

  # Other input checks.
  if (marker.callrate < 0 | marker.callrate > 1)
    stop("Specification of marker.callrate must be be between 0 and 1.")

  if (ind.callrate < 0 | ind.callrate > 1)
    stop("Specification of ind.callrate must be between 0 and 1.")

  if (maf < 0 | maf > 1)
    stop("Specification of maf must be between 0 and 1.")

  if (Fis < 0 | Fis > 1)
    stop("Specification of Fis must be between 0 and 1.")

  if (heterozygosity < 0 | heterozygosity > 1)
    stop("Specification of heterozygosity must be between 0 and 1.")

  # Check map if provided.
  if (!is.null(map)) {

    # Check map class.
    check.data_(data_ = "map", class_ = "data.frame")

    # Check map names.
    # Check mandatory variables in map.
    if(is.null(marker)){stop("\'marker' must be provided if \'map' is provided.")}

    # Check if they are present in the map.
    map.name.hit <- c(marker, chrom, pos) %in% names(map)
    if (!all(map.name.hit)){
      stop("Value provided to argument \'", c("marker", "chrom", "pos")[!map.name.hit],
           "' does not correspond to a variable in \'map'.")
    }

    # Match map and M.
    if (!identical(as.character(map[[marker]]), colnames(M))){
      stop("map[[marker]] and colnames(M) must be identical.")
    }
  }

  # # This is a slow test. Maybe not worth it. It is not necessary here.
  # if(!all(unique(c(M)) %in% c(0, 1, 2, na.string)) & message){
  #   message("Some of the values in M are not one of the following: ",
  #           paste0(c(0, 1, 2, na.string), collapse = ", "), ".")
  # }

  # Body ----------------------------------------------------------------------------------------

  # Initial info about the matrix.
  if (message) {
    message("Initial marker matrix M contains ", nrow(M),
            " individuals and ", ncol(M), " markers.")
  }

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

  # Check if all are compliant.
  if (!all(M %in% c(0, 1, 2, NA))) {
    stop("Data must be in numeric format: 0, 1, 2 and NA.")
  }

  # Remove markers with no valid information.
  miss.all <- colSums(is.na(M)) == nrow(M)

  if (any(miss.all)) {

    # Apply the removal.
    M <- M[, !miss.all]

    # Report.
    if (message){
      message("A total of ", sum(miss.all), " markers were removed for only having missing data.")
    }
  } ; rm(miss.all)


  # Generating some plots ---------------------------------------------------------------------------------

  if (plots){

    # Missing of individuals.
    # missingInd_DF <- data.frame(Ind =rowMeans(is.na(M)))
    missingInd_DF <- data.frame(Ind = (100 - callrate(M = M, margin = "row"))/100)
    missingInd_plot <- ggplot(missingInd_DF, aes(x=Ind)) +
      geom_histogram(fill='#0072B2', bins=40) +
      theme_classic() +
      xlab("Missing data per Individual")+
      ylab("Count")
    rm(missingInd_DF)

    # Missing of markers.
    # missingSNP_DF <- data.frame(SNP=colMeans(is.na(M)))
    missingSNP_DF <- data.frame(SNP = (100 - callrate(M = M, margin = "col"))/100)
    missingSNP_plot <- ggplot(missingSNP_DF, aes(x=SNP)) +
      geom_histogram(fill='#0072B2', bins=40) +
      theme_classic() +
      xlab("Missing data per SNP") +
      ylab("Count")

    # Histogram of MAF.
    qDF <- data.frame(MAF = maf(M = M))
    maf_plot <- ggplot(qDF, aes(x=MAF)) +
      geom_histogram(fill = '#0072B2', bins = 40) +
      theme_classic() +
      xlab("Minor Allele Frequency (MAF)")+
      ylab("Count")

    # Histogram of heterozygotes.
    het_DF <- data.frame(het = heterozygosity(M = M)[, "ho"])
    het_plot <- ggplot(het_DF, aes(x = het)) +
      geom_histogram(fill = '#0072B2', bins = 40) +
      theme_classic() +
      xlab("Heterozygotes")+
      ylab("Count")

    # Histogram of Fis.
    fis_DF <- data.frame(fis = abs(Fis(M = M)))
    fis_plot <- ggplot(fis_DF, aes(x = fis)) +
      geom_histogram(fill = '#0072B2', bins = 40) +
      theme_classic() +
      xlab("Fis")+
      ylab("Count")

  } else {
    missingInd_plot <- NULL
    missingSNP_plot <- NULL
    maf_plot <- NULL
    het_plot <- NULL
    fis_plot <- NULL
  }

  # Filtering markers -------------------------------------------------------------------------------------

  # Filtering process - 2 rounds (initializing objects).
  cr_mk_out <- 0 ; cr_id_out <- 0
  maf_out <- 0 ; fis_out <- 0
  h_out <- 0
  cr_mk_filter <- TRUE ; cr_id_filter <- TRUE
  maf_filter <- TRUE ; fis_filter <- TRUE
  h_filter <- TRUE

  for (blank_ in 1:2) {

    # Filtering markers by CR.
    if (marker.callrate < 1){
      cr_mk_filter <- 1 - callrate(M = M, margin = "col")/100 <= marker.callrate
      cr_mk_out <- cr_mk_out + sum(!cr_mk_filter)
      M <- M[, cr_mk_filter, drop = FALSE]
      rm(cr_mk_filter)
    }

    # Filtering callrate of individuals.
    if (ind.callrate < 1){
      cr_id_filter <- 1 - callrate(M = M, margin = "row")/100 <= marker.callrate
      cr_id_filter <- rowMeans(is.na(M)) <= ind.callrate
      cr_id_out <- cr_id_out + sum(!(cr_id_filter))
      M <- M[cr_id_filter, , drop = FALSE]
      rm(cr_id_filter)
    }

    # Filtering markers by MAF.
    if (maf > 0){
      q <- maf(M = M)
      maf_filter <- q > maf - .Machine$double.eps
      maf_out <- maf_out + sum(!maf_filter, na.rm = TRUE)
      M <- M[, maf_filter, drop = FALSE]
      rm(maf_filter, q)
    }

    # Filtering markers by heterozygosity.
    if (heterozygosity < 1){

      # Get observed heterozygosity.
      h <- heterozygosity(M = M)[, "ho"]

      # Get incidence vector.
      h_filter <- h <= heterozygosity

      # Add current run to the sum.
      h_out <- h_out + sum(!h_filter, na.rm = TRUE)

      # Apply filter.
      M <- M[, h_filter, drop = FALSE]

      # Remove objects.
      rm(h_filter, h)
    }

    if (Fis < 1){

      # Get Fis.
      fis <- Fis(M = M, margin = "col")

      # Get incidence vector.
      fis_filter <- abs(fis) <= Fis

      # Add current run to the sum.
      fis_out <- fis_out + sum(!fis_filter, na.rm = TRUE)

      # Apply filter.
      M <- M[, fis_filter, drop = FALSE]

      # Remove objects.
      rm(fis_filter, fis)
    }
  }

  # Some intermediate reporting.
  if (message){
    message("A total of ", cr_mk_out, " markers were removed because ",
            "their proportion of missing values was equal or larger than ",
            marker.callrate, ".")

    message("A total of ", cr_id_out, " individuals were removed because ",
            "their proportion of missing values was equal or larger than ",
            ind.callrate, ".")

    message("A total of ", maf_out, " markers were removed because ",
                 "their MAF was smaller than ", maf, ".")

    message("A total of ", h_out, " markers were removed because ",
                 "their heterozygosity was larger than ", heterozygosity, ".")

    message("A total of ", fis_out, " markers were removed because ",
                 "their |F| was larger than ", Fis, ".")

    missing.SNP <- sum(is.na(M))
    prop.miss <- 100*missing.SNP/(ncol(M)*nrow(M))
    message("Final cleaned marker matrix M contains ", round(prop.miss,2),
            "% of missing SNPs.")

    message("Final cleaned marker matrix M contains ", nrow(M),
            " individuals and ", ncol(M), " markers.")
    }

  # Simple mean imputation.
  if (impute){
    missing.SNP <- sum(is.na(M))
    if (missing.SNP == 0 & isTRUE(message)) {
      message('No imputation was performed as there are no missing marker data.')
    } else {

      # Loop through markers and impute with mean.
      for(i in 1:ncol(M)){
        M[is.na(M[,i]), i] <- mean(M[,i], na.rm=TRUE)
      }

      # Polishing the dataset.
      M <- round(M, digits)
      if (isTRUE(message)){
        prop.miss <- 100*missing.SNP/(ncol(M)*nrow(M))
        message("A total of ", missing.SNP, " missing values were imputed, ",
                "corresponding to ", round(prop.miss,2), "% of the total number of SNPs.")
      }
    }
  }

  # Finalize ------------------------------------------------------------------------------------

  # Remove eventual cleaned markers from M in the map.
  if (!is.null(map)){

    # Get match index.
    matches <- na.omit(match(colnames(M), as.character(map[[marker]])))

    # Applies to map; This separation is done because of data.table.
    map <- map[matches,]
  }

  return(list(M.clean = M, map = map, plot.missing.ind = missingInd_plot,
         plot.missing.SNP = missingSNP_plot, plot.heteroz = het_plot, plot.Fis = fis_plot,
         plot.maf = maf_plot))

}
