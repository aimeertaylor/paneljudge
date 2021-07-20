###########################################################################
#' Function to estimate relatedness and switch rate parameters
#'
#' Given a matrix of marker allele frequencies, a vector of inter-marker
#' distances, and a matrix of genotype calls for a pair of haploid genotypes,
#' \code{estimate_r_and_k} returns the maximum likelihood estimates of the
#' relatedness parameter, \eqn{r}, and the switch rate parameter, \eqn{k}, under
#' the HMM described in [1].
#'
#' @param fs Matrix of marker allele frequencies, i.e. the \eqn{ft}s in [1].
#'   Specifically, a \eqn{m} by \eqn{Kmax} matrix, where \eqn{m} is the marker
#'   count and \eqn{Kmax} is the maximum cardinality (per-marker allele count)
#'   observed over all \eqn{m} markers. If, for any \eqn{t = 1,...,m}, the
#'   maximum cardinality exceeds that of the \eqn{t}-th marker (i.e. if
#'   \eqn{Kmax > Kt}), then all \code{fs[t,1:Kt]} are in (0,1] and all
#'   \code{fs[t,(Kt+1):Kmax]} are zero. For example, if \eqn{Kt = 2} and
#'   \eqn{Kmax = 4} then \code{fs[t,]} might look like \code{[0.3, 0.7, 0, 0]}.
#' @param ds Vector of \eqn{m} inter-marker distances, i.e. the \eqn{dt}s in
#'   [1]. The \eqn{t}-th element of the inter-marker distance vector,
#'   \code{ds[t]}, contains the distance between marker \eqn{t} and \eqn{t+1}
#'   such that \code{ds[m] = Inf}, where \eqn{m} is the marker count. (Note that
#'   this differs slightly from [1], where \code{ds[t]} contains the distance
#'   between marker \eqn{t-1} and \eqn{t}). Distances between markers on
#'   different markers are also considered infinite, i.e. if the chromosome of
#'   marker \eqn{t+1} is not equal to the chromosome of \eqn{t}-th marker,
#'   \code{ds[t] = Inf}.
#' @param Ys Matrix of genotypes calls for a pair of simulated haploid
#'   genotypes, i.e. the \eqn{Yt}s of the \eqn{i}-th and \eqn{j}-th haploid
#'   genotypes in [1]. Specifically, a \eqn{m} by 2 matrix, where \eqn{m} is the
#'   marker count and each column contains a haploid genotype. For all \eqn{t =
#'   1,...,m} markers, alleles are enumerated 0 to \eqn{Kt-1}, where \eqn{Kt} is
#'   the cardinality (per-marker allele count) of the \eqn{t}-th marker. For
#'   example, if \eqn{Kt = 2}, both \code{Ys[t,1]} and \code{Ys[t,2]} are either
#'   0 or 1.
#' @param epsilon Genotyping error, i.e. \eqn{\epsilon} in [1]. The genotyping
#'   error is the probability of miscalling one specific allele for another. As
#'   such, the error rate for the t-th marker, \eqn{(Kt-1)\epsilon}, scales with
#'   \eqn{Kt} (the per-marker allele count, cardinality).
#' @param rho Recombination rate, i.e. \eqn{\rho} in [1]. The recombination rate
#'   corresponds to the probability of a crossover per base pair. It is assumed
#'   constant across the genome under the HMM of [1]. Its default value
#'   corresponds to an average rate estimated for \emph{Plasmodium falciparum}
#'   [2].
#' @param kinit Switch rate parameter value used to initialise optimization of
#'   the negative loglikelihood.
#' @param rinit Relatedness parameter value used to initialise optimization of
#'   the negative loglikelihood.
#' @param warn_fs Logical indicating if the function should return warnings
#'   following allele frequency checks.
#'
#' @return Maximum likelihood estimates of the switch rate parameter, \eqn{k},
#'   and relatedness parameter, \eqn{r}.
#'
#' @examples
#' # First stimulate some data
#' simulated_Ys <- simulate_Ys(fs = frequencies$Colombia, ds = markers$distances, k = 5, r = 0.25)
#'
#' # Second estimate the switch rate parameter, k, and relatedness parameter, r
#' estimate_r_and_k(fs = frequencies$Colombia, ds = markers$distances, Ys = simulated_Ys)
#'
#' @references \enumerate{ \item Taylor, A.R., Jacob, P.E., Neafsey, D.E. and
#'   Buckee, C.O., 2019. Estimating relatedness between malaria parasites.
#'   Genetics, 212(4), pp.1337-1351. \item Miles, A., Iqbal, Z., Vauterin, P.,
#'   Pearson, R., Campino, S., Theron, M., Gould, K., Mead, D., Drury, E.,
#'   O'Brien, J. and Rubio, V.R., 2016. Indels, structural variation, and
#'   recombination drive genomic diversity in Plasmodium falciparum. Genome
#'   research, 26(9), pp.1288-1299.}
#'
#' @export
###########################################################################
estimate_r_and_k <- function(fs, ds, Ys, epsilon = 0.001, rho = 7.4 * 10 ^ (-7),
                             kinit = 50, rinit = 0.5, warn_fs = TRUE) {

  # Convert to a into matrix if not already (loglikelihood_cpp expects a matrix)
  # and perform some checks
  if (!is.matrix(fs)) fs <- as.matrix(fs)
  # If no errors, returns numeric logic for fs zero/non-zero
  fs_checks_return <- fs_checks(fs, warn = warn_fs, do_return = TRUE)

  # Check Ys for NAs
  if(any(is.na(Ys))) stop("Missing values detected in Ys.\n  Please remove and recompute ds accordingly.")

  # Check for alleles that are permissible only if epsilon exceeds zero.
  # These alleles will break the hmmloglikelihood.cpp code if epsilon is zero.
  problem <- "Some per-marker allele counts exceed per-marker non-zero allele frequencies."
  fs_checks_return$fs01 <- 1 * (fs > fs_checks_return$non_zero_fs_lb)
  Kts <- rowSums(fs_checks_return$fs01)
  if (any(sapply(1:nrow(Ys), function(i) Ys[i] > (Kts[i]-1)))) {
    if (epsilon > 0) warning (paste0(problem, " Data are permissible due to non-zero epsilon."))
    else if (epsilon == 0)  stop(paste0(problem, " Data are incompatible with zero-valued epsilon."))
  }

  # Define the function to pass to optim()
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, fs, ds, epsilon, rho)

  # Optimise the negative log likelihood
  optimization <- optim(par = c(kinit, rinit), fn = function(x) - ll(x[1], x[2]))

  # Extract and name estimates
  rkhats <- c("khat" = optimization$par[1], "rhat" = optimization$par[2])

  if (all(rkhats == c(kinit, rinit))) {
    warning("optimization has returned initial parameter values. Data are possibly uniformative.")
  }

  # End of function
  return(rkhats)
}
