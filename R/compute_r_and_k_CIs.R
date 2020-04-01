#' @importFrom doRNG %dorng%
#' @examples
#' compute_r_and_k_CIs(frequencies = frequencies$Colombia, distances = markers$dt, khat, rhat)

compute_r_and_k_CIs <- function(frequencies, distances, khat, rhat,
                                confidence = 95, nboot = 100, do_not_use_core_count = 2,
                                epsilon = 0.001, rho = 7.4 * 10^(-7), kinit = 50, rinit = 0.5){

  # Register available cores in order to run code in parallel
  # Note that packages parallel and foreach are available through doParallel, which is imported
  doParallel::registerDoParallel(cores = parallel::detectCores() - do_not_use_core_count)

  # Generate parametric boot strap r and k estimates
  # %dorng% is similar to %dopar%, however
  # %dorng% loops are reproducible whereas %dopar% loops are not,
  # %dorng% returns the whole sequence of RNG seeds as an attribute whereas %dopar% does not.
  attr(r1, 'rng')
  rkhats_boot = foreach::foreach(iboot = 1:nboot, .combine = rbind) %dorng% {
    simulated_genotype_pair <- simulate_data(frequencies, distances, k = khat, r = rhat, epsilon, rho)
    rkhats_boot <- estimate_r_and_k(frequencies, distances, Ys = simulated_genotype_pair, epsilon, rho, kinit, rinit)
  }

  # Calculata CIs as quantiles
  alpha <- (1-confidence/100)
  CI <- apply(rkhats_boot, 2, quantile, probs = c(alpha/2, 1-alpha/2))
  if(any(is.na(CI))){stop('Bootstrap CI is NA')}

  # Return CIs
  return(CI)
}



