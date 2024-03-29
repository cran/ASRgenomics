% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synthetic_cross.R
\name{synthetic.cross}
\alias{synthetic.cross}
\title{Generates a molecular matrix M for hypothetical crosses based on the
genomic information of the parents}
\usage{
synthetic.cross(
  M = NULL,
  ped = NULL,
  indiv = NULL,
  mother = NULL,
  father = NULL,
  heterozygote.action = c("useNA", "exact", "fail", "expected"),
  na.action = c("useNA", "expected"),
  message = TRUE
)
}
\arguments{
\item{M}{A matrix with marker data of full form (\eqn{n \times p}), with \eqn{n} individuals
(mothers and fathers) and \eqn{p} markers.
Individual and marker names are assigned to \code{rownames} and \code{colnames}, respectively.
Data in matrix is coded as 0, 1, 2 (integer or numeric) (default = \code{NULL}).}

\item{ped}{A data frame with three columns containing only the pedigree of the hypothetical offspring.
(not pedigree of parents)
It should include the three columns for individual, mother and father (default = \code{NULL}).}

\item{indiv}{A character indicating the column in \code{ped} data frame containing the identification
of the offspring (default = \code{NULL}).}

\item{mother}{A character indicating the column in \code{ped} data frame containing the identification
of the mother (default = \code{NULL}).}

\item{father}{A character indicating the column in \code{ped} data frame containing the identification
of the father (default = \code{NULL}).}

\item{heterozygote.action}{Indicates the action to take when heterozygotes are found in a marker.
Options are: \code{"useNA"}, \code{"exact"}, \code{"fail"}, and \code{"expected"}.
See details for more information (default = \code{"useNA"})}

\item{na.action}{Indicates the action to take when missing values are found in a marker.
Options are: \code{"useNA"} and \code{"expected"}.
See details for more information (default = \code{"useNA"}).}

\item{message}{If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).}
}
\value{
A molecular matrix \eqn{\boldsymbol{M}} containing the genotypes generated/imputed for the
hypothetical cross.
}
\description{
This function generates (or imputes) a molecular matrix for offspring
from hypothetical crosses based on the genomic information from the parents.
This is a common procedure in species such as maize, where only the parents
(inbred lines) are genotyped, and this information is used to generate/impute
the genotypic data of each of the hybrid offspring.
This function can be also used for bulked DNA analyses, in order to obtain an
bulked molecular matrix for full-sib individuals were only parents are genotyped.
}
\details{
For double-haploids, almost the totality of the markers (except for genotyping errors)
will be homozygotic reads. But in many other cases (including recombinant inbred lines)
there will be a proportion of heterozygotic reads. In these case, it is
very difficult to infer (impute) the exact genotype of a given offspring individual.
For example, if parents are 0 (AA) and 1 (AC) then offsprings will differ given this
Mendelian sampling. However, different strategies exist to determine the
expected value for that specific cross (if required), which are detailed below
using the option \code{heterozygote.action}.

\itemize{

\item{If \code{heterozygote.action =  "useNA"},
the generated offspring will have, for the heterozygote read, an \code{NA},
and no markers are removed.
Hence, no attempt will be done to impute/estimate this value.
}

\item{If \code{heterozygote.action = "exact"},
any marker containing one or more heterozygote reads will be removed.
Hence, inconsistent markers are fully removed from the \eqn{\boldsymbol{M}} matrix.
}

\item{If \code{heterozygote.action = "fail"},
function stops and informs of the presence of heterozygote reads.
}

\item{If \code{heterozygote.action = "expected"},
then an algorithm is implemented, on the heterozygote read to determine its
expected value. For example, if parents are 0 and 1, then the expected value
(with equal probability) is 0.5. For a cross between two heterozygotes,
the expected value is: \eqn{0(1/4) + 1(1/2) + 2(1/4) = 1}. And for a cross
between 1 and 2, the expected value is: \eqn{1(1/2) + 2(1/2) = 1.5}}
}

Missing value require special treatment, and an imputation strategy is detailed
below as indicated using the option \code{na.action}.

\itemize{

\item{If \code{na.action = "useNA"}, if at least one of the parental reads
is missing values for a given marker then it will be assigned as missing for
the hypothetical cross. Hence, no attempt will be done to impute/estimate
this value.}

\item{If \code{na.action = "expected"}, then an algorithm is implemented that
will impute the expected read of the cross if the genotype of \strong{one of
the parents is missing} (\emph{e.g.}, cross between 0 and NA). Calculations
are based on parental allelic frequencies \eqn{p} and \eqn{q} for the given
marker. The expressions for expected values are detailed below.}

\itemize{

\item{If the genotype of the non-missing parent read is 0.

\eqn{q^2} (probability that the missing parent is 0) x 0 (expected value of the offspring from a 0 x 0 cross: \eqn{0(1/1)}) +

\eqn{2pq} (probability that the missing parent is 1) x 0.5 (expected offspring from a 0 x 1 cross: \eqn{0(1/2) + 1(1/2)}) +

\eqn{q^2} (probability that the missing parent is 2) x 1 (expected offspring from a 0 x 2 cross: \eqn{1(1/1)})}

\item{If the genotype of the non-missing parent read is 1.

\eqn{q^2} (probability that the missing parent is 0) x 0.5 (offspring: \eqn{0(1/2) + 1(1/2)}) +

\eqn{2pq} (probability that the missing parent is 1) x 1 (offspring: \eqn{0(1/4) + 1(1/2) + 2(1/4)}) +

\eqn{q^2} (probability that the missing parent is 2) x 1.5 (offspring: \eqn{1(1/2) + 2(1/2)})}

\item{If the genotype of the non-missing parent read is 2.

\eqn{q^2} (probability that the missing parent is 0) x 1 (offspring: \eqn{1(1/1)}) +

\eqn{2pq} (probability that the missing parent is 1) x 1.5 (offspring: \eqn{1(1/2) + 2(1/2)}) +

\eqn{q^2} (probability that the missing parent is 2) x 2 (offspring: \eqn{2(1/1)})}
}

Similarly, the calculation of the expected read of a cross when \strong{both parents are missing} is
also based on population allelic frequencies for the given marker.
The expressions for expected values are detailed below.

              \eqn{q^2 \times q^2} (probability that both parents are 0) x 0 (expected value of the offspring from a 0 x 0 cross: 0(1/1)) +

              \eqn{2 \times (q^2 \times 2pq)} (probability that the first parent is 0 and the second is 1; this requires
the multiplication by 2 because it is also possible that the first parent is 1 and the second is 0)
x 0.5 (offspring: \eqn{0(1/2) + 1(1/2)}) +

              \eqn{2 \times (q^2 \times p^2)} (this could be 0 x 2 or 2 x 0) x 1 (offspring: \eqn{1(1/1)}) +

              \eqn{2pq \times 2pq} (both parents are 1) x 1 (offspring: \eqn{0(1/4) + 1(1/2) + 2(1/4)}) +

              \eqn{2 \times (2pq \times q2)} (this could be 1 x 2 or 2 x 1) x 1.5 (offspring: \eqn{1(1/2) + 2(1/2)}) +

              \eqn{p^2 \times p^2} (both parents are 2) x 2 (offspring: \eqn{2(1/1)})

Note that the use of \code{na.action = "expected"} is recommended when
a large number of offspring will conform the hybrid cross (such as
with bulked DNA analyses) for family groups with reasonable number of individuals.

\strong{Warning}. If \code{"expected"} is used for \code{heterozygote.action} or \code{na.action},
direct transformation of the molecular data to other codings (\emph{e.g.},
dominance matrix coded as \code{c(0,1,0)}) is not recommended.
}
}
\examples{
# Create dummy pedigree (using first 10 as parents).
ped <- data.frame(
 male = rownames(geno.apple)[1:5],
 female = rownames(geno.apple)[6:10])
ped$offs <- paste(ped$male, ped$female, sep = "_")
ped

# Select portion of M for parents.
Mp <- geno.apple[c(ped$male, ped$female), 1:15]

# Get genotype of crosses removing markers with heterozygotes.
synthetic.cross(
 M = Mp, ped = ped,
 indiv = "offs", mother = "female", father = "male",
 heterozygote.action = "exact",
 na.action = "useNA")

# Request the synthetic cross to be NA in the respective samples.
synthetic.cross(
 M = Mp, ped = ped,
 indiv = "offs", mother = "female", father = "male",
 heterozygote.action = "useNA",
 na.action = "useNA")

# Get genotype of crosses and use expected values.
synthetic.cross(
 M = Mp, ped = ped,
 indiv = "offs", mother = "female", father = "male",
 heterozygote.action = "expected", na.action = "expected")

}
