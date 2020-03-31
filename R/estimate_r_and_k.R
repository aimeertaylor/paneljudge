## Mechanism to compute MLEs given fs, Ys, epsilon and distances
estimate_r_and_k <- function(frequencies, distances, Ys,
                             epsilon = 0.001, , rho = 7.4 * 10^(-7), kinit = 50, rinit = 0.5){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho)
  optimization <- optim(par = c(kinit, rinit), fn = function(x) - ll(x[1], x[2]))
  rkhats <- optimization$par
  return(rkhats)
}
