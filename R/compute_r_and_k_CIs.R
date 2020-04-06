###########################################################################
#' Compute confidence intervals for relatedness and switch rate parameters
#'
#' For specified confidence (default 95\%), given a matrix of input frequencies
#' and a vector of marker distances, a switch rate parameter estimate, khat, and
#' a relatedness parameter estimate, rhat. \code{compute_r_and_k_CIs()} returns
#' the confidence intervals around the parameter estimates. Confidence intervals
#' are approximate. They are approximated using the parametric bootstrap draws
#' under the HMM described in Taylor et al. 2009 The parametric bootstrap draws
#' are generated in parallel using a specifying number of cores. The quality and
#' compute time of the approximation scales with nboot.
#'
#' @param frequencies ndata (marker count) by Kmax (max cardinality over m
#'   markers) matrix if Kt < Kmax for any t in 1:m, then fs[t,1:Kt] in (0,1) &
#'   fs[t,(Kt+1):Kmax] = 0 Example: if Kt = 2 < Kmax = 4 then fs[t,] might look
#'   like [0.2, 0.7, 0, 0].
#' @param distances vector where distances[t] contains the distance between
#'   position t and t+1 or equivalently, gendist[t-1] contains the distance
#'   between position t-1 and t, for p > 1. s.t. only m-1 first entries of
#'   gendist are being used.
#' @param confidence confidence level (percentage) of confidence interval
#'   (default 95\%)
#' @param nboot number of parametric bootstrap estimates from which to compute
#'   the CI. Larger values provide a better approximation but CIs will take
#'   longer to compute.
#' @param do_not_use_core_count number of cores to keep free on computer
#'   (defaults to two)
#' @param epsilon genotyping error (probability of miscall).
#' @param rho recombination rate in probability of break point per base pair.
#' @param kinit switch rate parameter value used to initialise loglikelihood
#'   optimization.
#' @param rinit relatedness parameter value used to initialise loglikelihood
#'   optimization.
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#' @return Confidence interval bounds around input switch rate parameter, k, and
#'   relatedness parameter, r.
#' @examples
#' # First, stimulate some data
#' simulated_Ys <- simulate_Ys(fs = frequencies$Colombia, ds = markers$dt, k = 5, r = 0.25)
#'
#' # Second, estimate the switch rate parameter, k, and relatedness parameter, r
#' krhat <- estimate_r_and_k(fs = frequencies$Colombia, ds = markers$dt, Ys = simulated_Ys)
#'
#' # Third, compute confidence intervals (CIs)
#' compute_r_and_k_CIs(frequencies$Colombia, markers$dt, khat = krhat['khat'], rhat = krhat['rhat'])
#'
#' @references Taylor, A.R., Jacob, P.E., Neafsey, D.E. and Buckee, C.O., 2019.
#'   Estimating relatedness between malaria parasites. Genetics, 212(4),
#'   pp.1337-1351.
#' @export
###########################################################################

compute_r_and_k_CIs <- function(frequencies, distances, khat, rhat,
                                confidence = 95, nboot = 100, do_not_use_core_count = 2, ...){

  # Retrieve all additional parameters
  all_params <- list(...)

  # Define a subset of additional parameters to pass to simulate_Ys()
  # Alternatively, could use R.utils::doCall() with .ignoreUnusedArgs=TRUE
  # instead of base::do.call(); see call to simulate_Ys() below.
  sim_params <- all_params
  sim_params$kinit <- NULL
  sim_params$rinit <- NULL

  # Register available cores in order to run code in parallel
  # Note that packages parallel and foreach are available through doParallel, which is imported
  doParallel::registerDoParallel(cores = parallel::detectCores() - do_not_use_core_count)

  # Generate parametric boot strap r and k estimates
  # %dorng% is similar to %dopar%, however
  # %dorng% loops are reproducible whereas %dopar% loops are not,
  # %dorng% returns the whole sequence of RNG seeds as an attribute whereas %dopar% does not.
  rkhats_boot = foreach::foreach(iboot = 1:nboot, .combine = rbind) %dorng% {
    simulated_genotype_pair <- do.call(simulate_Ys, args = c(list(frequencies, distances, k = khat, r = rhat), sim_params))
    rkhats_boot <- estimate_r_and_k(frequencies, distances, Ys = simulated_genotype_pair, ...)
  }

  # Calculata CIs as quantiles
  alpha <- (1-confidence/100)
  CI <- apply(rkhats_boot, 2, quantile, probs = c(alpha/2, 1-alpha/2))
  if(any(is.na(CI))){stop('Bootstrap CI is NA')}

  # Return CIs
  return(t(CI))
}



