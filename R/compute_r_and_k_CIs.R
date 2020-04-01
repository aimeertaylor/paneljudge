rhats_hmm_boot = as.numeric(foreach(iboot = 1:nboot, .combine = c) %dorng% {
  Ys_boot <- simulate_Ys_hmm(frequencies, distances = dataset_$dt, k = mle_df$khat_hmm[ind], r = mle_df$rhat_hmm[ind], epsilon)
  ll_hmm <- function(k, r) loglikelihood_cpp(k, r, Ys_boot, frequencies, dataset_$dt, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optim(par = c(50, 0.5), fn = function(x) - ll_hmm(x[1], x[2]))
  optimization$par[2]
})

CI <- as.numeric(quantile(rhats_hmm_boot, probs = c(0.025, 0.975)))
if(any(is.na(CI))){stop('Bootstrap CI is NA')}
