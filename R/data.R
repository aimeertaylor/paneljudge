#' Markers of the example GTseq panel.
#'
#' A data set of marker attributes for markers pertaining to the example GTseq panel.
#'
#' @format A data frame with 126 rows and 8 variables:
#' \describe{
#'   \item{Amplicon_name}{Name (character) of the microhaplotype marker ("Amplicon" because typed using an amplicon)}
#'   \item{Chr}{Chromosome (character) of the microhaplotype marker}
#'   \item{Start}{First base pair (integer) of the microhaplotype marker}
#'   \item{Stop}{Last base pair  (integer) of the microhaplotype marker}
#'   \item{length}{Length (integer) of the microhaplotype marker in base pairs}
#'   \item{pos}{Mid-point (numeric) of the microhaplotype marker}
#'   \item{chrom}{Chromosome (numeric) of the microhaplotype marker}
#'   \item{distance}{Inter mid-point distance (numeric) between the microhaplotype marker and its subsequent neighbour}
#' }
#' @source see \url{https://github.com/artaylor85/paneljudge/blob/master/data_raw/Process_GTseq.R}
"markers"

#' Allele frequencies of the example GTseq panel.
#'
#' A data set of allele frequencies for four countries: Colombia, French Guiana, Mali and Sengal.
#'
#' @format Each entry of the list is a matrix, \code{fs} say, with \eqn{m=126}
#'   rows and \eqn{Kmax=44} variables, where \eqn{m} is the marker count and
#'   \eqn{Kmax} is the maximum cardinality (per-marker allele count) observed
#'   over all \eqn{m} markers. If, for any \eqn{t = 1,...,m}, the maximum
#'   cardinality exceeds that of the \eqn{t}-th marker (i.e. if \eqn{Kmax >
#'   Kt}), then all \code{fs[t,1:Kt]} are in (0,1) and all
#'   \code{fs[t,(Kt+1):Kmax]} are zero. For example, for PF3D7_0103600 in Colombia, \eqn{Kt
#'   = 2} and \code{frequencies$Colombia["PF3D7_0103600",] = (0.687075, 0.312925, 0, ..., 0)}.
#'   \describe{
#'   \item{Allele.1}{Frequency (numeric) of the first allele}
#'   ...
#'   \item{Allele.44}{Frequency (numeric) of the \eqn{Kmax} allele}
#' }
#' @source see \url{https://github.com/artaylor85/paneljudge/blob/master/data_raw/Process_GTseq.R}
"frequencies"

#' Chromosome lengths.
#'
#' Lengths in base pairs of chromosomes Pf3D7_01_v3 to Pf3D7_14_v3 of the 3D7
#' Plasmodium falciparum reference genome listed on PlasmoDB (see url below).
#'
#' @format A numeric vector named by the chromosome number.
#'
#' @source \url{https://plasmodb.org/plasmo/showApplication.do}
"chr_lengths"

