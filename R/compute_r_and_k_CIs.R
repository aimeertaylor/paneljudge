###########################################################################
#' Compute confidence interval bounds for relatedness and switch rate parameters
#'
#' Given a matrix of marker allele frequencies, a vector of inter-marker
#' distances, and estimates of the relatedness and switch rate parameters,
#' \code{compute_r_and_k_CIs()} returns confidence interval bounds around the
#' parameter estimates. The default confidence is 95\%. The interval bounds are
#' approximate. They are generated using parametric bootstrap draws of the
#' parameter estimates based on genotype calls for haploid genotype pairs
#' simulated under the HMM described in [1] using the input parameter estimates.
#' The quality of the approximation increases and compute time scales with the
#' number of parametric bootstrap draws, which are generated in parallel using a
#' specified number of cores.
#'
#' @param fs Matrix of marker allele frequencies, i.e. the \eqn{ft}s in [1].
#'   Specifically, a \eqn{m} by \eqn{Kmax} matrix, where \eqn{m} is the marker
#'   count and \eqn{Kmax} is the maximum cardinality (per-marker allele count)
#'   observed over all \eqn{m} markers. If, for any \eqn{t = 1,...,m}, the
#'   maximum cardinality exceeds that of the \eqn{t}-th marker (i.e. if
#'   \eqn{Kmax > Kt}), then all \code{fs[t,1:Kt]} are in (0,1) and all
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
#' @param khat Estimate of the switch rate parameter, i.e. estimate of \eqn{k} in [1].
#' @param rhat Estimate of the relatedness parameter, i.e. estimate of \eqn{r} in [1].
#' @param confidence Confidence level (percentage) of confidence interval
#'   (default 95\%).
#' @param nboot Number of parametric bootstrap draws from which to compute
#'   the confidence interval bounds. Larger values provide a better approximation
#'   but prolong computation.
#' @param core_count Number of cores to use to do computation.
#' Set to 2 or more for parallel computation.
#' Defaults to the number detected on the machine minus one.
#' @param ... Arguments to be passed to \code{simulate_Ys()} and \code{estimate_r_and_k()}.
#'
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#'
#' @return Confidence interval bounds around input switch rate parameter, \eqn{k}, and
#'   relatedness parameter, \eqn{r}.
#'
#' @examples
#' # First, stimulate some data
#' simulated_Ys <- simulate_Ys(fs = frequencies$Colombia, ds = markers$distances, k = 5, r = 0.25)
#'
#' # Second, estimate the switch rate parameter, k, and relatedness parameter, r
#' krhat <- estimate_r_and_k(fs = frequencies$Colombia, ds = markers$distances, Ys = simulated_Ys)
#'
#' # Third, compute confidence intervals (CIs)
#' compute_r_and_k_CIs(frequencies$Colombia, markers$distances, khat = krhat['khat'], rhat = krhat['rhat'])
#'
#' @references \enumerate{ \item Taylor, A.R., Jacob, P.E., Neafsey, D.E. and Buckee, C.O., 2019.
#'   Estimating relatedness between malaria parasites. Genetics, 212(4),
#'   pp.1337-1351.}
#'
#' @export
###########################################################################

compute_r_and_k_CIs <- function(fs, ds, khat, rhat, confidence = 95, nboot = 100,
                                core_count = parallel::detectCores()-1, ...){

  # Retrieve all additional parameters
  all_params <- list(...)

  # Define a subset of additional parameters to pass to simulate_Ys() using base::do.call(); see below.
  # Alternatively, could use R.utils::doCall() with .ignoreUnusedArgs=TRUE.
  sim_params <- all_params
  sim_params$kinit <- NULL
  sim_params$rinit <- NULL

  # Register available cores in order to run code in parallel
  # Note that packages parallel and foreach are available through doParallel, which is imported
  doParallel::registerDoParallel(cores = core_count)

  # Generate parametric bootstrap draw of r and k estimates
  # %dorng% is similar to %dopar%, however
  # %dorng% loops are reproducible whereas %dopar% loops are not,
  # %dorng% returns the whole sequence of RNG seeds as an attribute whereas %dopar% does not.
  rkhats_boot = foreach::foreach(iboot = 1:nboot, .combine = rbind) %dorng% {
    simulated_Ys <- do.call(simulate_Ys, args = c(list(fs, ds, k = khat, r = rhat), sim_params))
    rkhats_boot <- estimate_r_and_k(fs, ds, Ys = simulated_Ys, ...)
  }

  # Calculate CIs as quantiles of the parametric bootstrap draw
  alpha <- (1-confidence/100)
  CI <- apply(rkhats_boot, 2, quantile, probs = c(alpha/2, 1-alpha/2))
  if (any(is.na(CI))) stop('Bootstrap CI is NA')

  # End of function
  return(t(CI))
}



