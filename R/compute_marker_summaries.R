###########################################################################
#' Function to compute marker effective cardinalities
#'
#' Given a matrix of marker allele frequencies,
#' \code{compute_eff_cardinalities} returns the effective cardinalites of
#' \eqn{t = 1,...,m} markers, where \eqn{m} is the marker count. Effective
#' cardinalities are per-marker allele counts that account for inequifrequenct
#' alleles. Each effective cardinality is calculated as described in [1], i.e.
#' without correcting for finite sample sizes or considering uncertainty.
#'
#' @param fs Matrix of marker allele frequencies, i.e. the \eqn{ft}s in [1].
#'   Specifically, a \eqn{m} by \eqn{Kmax} matrix, where \eqn{m} is the marker
#'   count and \eqn{Kmax} is the maximum cardinality (per-marker allele count)
#'   observed over all \eqn{m} markers. If, for any \eqn{t = 1,...,m}, the
#'   maximum cardinality exceeds that of the \eqn{t}-th marker (i.e. if
#'   \eqn{Kmax > Kt}), then all \code{fs[t,1:Kt]} are in (0,1] and all
#'   \code{fs[t,(Kt+1):Kmax]} are zero. For example, if \eqn{Kt = 2} and
#'   \eqn{Kmax = 4} then \code{fs[t,]} might look like \code{[0.3, 0.7, 0, 0]}.
#' @param warn_fs Logical indicating if the function should return warnings
#'   following allele frequency checks.
#'
#' @return Effective cardinalities for \eqn{t = 1,\ldots,m} markers.
#'
#' @examples
#' compute_eff_cardinalities(fs = frequencies$Colombia)
#'
#' @references \enumerate{ \item Taylor, A.R., Jacob, P.E., Neafsey, D.E. and
#'   Buckee, C.O., 2019. Estimating relatedness between malaria parasites.
#'   Genetics, 212(4), pp.1337-1351.}
#'
#' @export
###########################################################################
compute_eff_cardinalities <- function(fs, warn_fs = TRUE) {
  fs_checks(fs, warn = warn_fs) # Check frequencies
  eff_cardinalities <- 1 / rowSums(fs^2)
  return(eff_cardinalities)
}


###########################################################################
#' Function to compute marker diversities
#'
#' Given a matrix of marker allele frequencies, \code{compute_diversities}
#' returns the diversities of \eqn{t = 1,...,m} markers, where \eqn{m} is the
#' marker count. Each diversity is calculated as described in [1], i.e. without
#' correcting for finite sample sizes or considering uncertainty.
#'
#' @param fs Matrix of marker allele frequencies, i.e. the \eqn{ft}s in [1].
#'   Specifically, a \eqn{m} by \eqn{Kmax} matrix, where \eqn{m} is the marker
#'   count and \eqn{Kmax} is the maximum cardinality (per-marker allele count)
#'   observed over all \eqn{m} markers. If, for any \eqn{t = 1,...,m}, the
#'   maximum cardinality exceeds that of the \eqn{t}-th marker (i.e. if
#'   \eqn{Kmax > Kt}), then all \code{fs[t,1:Kt]} are in (0,1) and all
#'   \code{fs[t,(Kt+1):Kmax]} are zero. For example, if \eqn{Kt = 2} and
#'   \eqn{Kmax = 4} then \code{fs[t,]} might look like \code{[0.3, 0.7, 0, 0]}.
#' @param warn_fs Logical indicating if the function should return warnings
#'   following allele frequency checks.
#'
#' @return Diversities for \eqn{t = 1,\ldots,m} markers.
#'
#' @examples
#' compute_diversities(fs = frequencies$Colombia)
#'
#' @references \enumerate{ \item Taylor, A.R., Jacob, P.E., Neafsey, D.E. and
#'   Buckee, C.O., 2019. Estimating relatedness between malaria parasites.
#'   Genetics, 212(4), pp.1337-1351.}
#'
#' @export
###########################################################################
compute_diversities <- function(fs, warn_fs = TRUE) {
  fs_checks(fs, warn = warn_fs) # Check frequencies
  diversities <- 1 - rowSums(fs^2)
  return(diversities)
}
