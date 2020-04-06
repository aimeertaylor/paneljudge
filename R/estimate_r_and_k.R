###########################################################################
#' Estimate relatedness and switch rate parameters
#'
#' Given a matrix of frequencies fs and a vector of marker ds,
#' \code{estimate_r_and_k()} returns the maximum likelihood estimates of the
#' relatedness parameter, r, and the switch rate parameter, k, under the HMM
#' described in [1].
#'
#' @param fs ndata (marker count) by Kmax (max cardinality over m
#'   markers) matrix if Kt < Kmax for any t in 1:m, then fs[t,1:Kt] in (0,1) &
#'   fs[t,(Kt+1):Kmax] = 0 Example: if Kt = 2 < Kmax = 4 then fs[t,] might look
#'   like [0.2, 0.7, 0, 0].
#' @param ds vector where ds[t] contains the distance
#'   between position t and t+1 or equivalently, gendist[t-1] contains the
#'   distance between position t-1 and t, for p > 1. s.t. only m-1 first entries
#'   of gendist are being used.
#' @param Ys Data on a pair of haploid genotypes
#' @param epsilon genotyping error (probability of miscall).
#' @param rho recombination rate in probability of break point per base pair.
#' @param kinit switch rate parameter value used to initialise loglikelihood
#'   optimization.
#' @param rinit relatedness parameter value used to initialise loglikelihood
#'   optimization.
#'
#' @return Maximum likelihood estimates (MLEs) of the switch rate parameter, k, and relatedness parameter, r.
#' @examples
#' # First stimulate some data
#' Simulated_genotype_pair <- simulate_data(fs = frequencies$Colombia, gendist = markers$dt, k = 5, r = 0.25)
#' # Second estimate the switch rate parameter, k, and relatedness parameter, r
#' estimate_r_and_k(fs = frequencies$Colombia, ds = markers$dt, Ys = Simulated_genotype_pair)
#'
#' @references Taylor, A.R., Jacob, P.E., Neafsey, D.E. and Buckee, C.O., 2019.
#'   Estimating relatedness between malaria parasites. Genetics, 212(4),
#'   pp.1337-1351.
#' @export
###########################################################################
estimate_r_and_k <- function(fs, ds, Ys, epsilon = 0.001, rho = 7.4 * 10^(-7), kinit = 50, rinit = 0.5){

  # Convert to a into matrix if not already (loglikelihood_cpp expects a matrix)
  if (!is.matrix(fs)) {fs <- as.matrix(fs)}

  # Define the function to pass to optim
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, fs, ds, epsilon, rho)

  # Optimise the negative log likelihood
  optimization <- optim(par = c(kinit, rinit), fn = function(x) - ll(x[1], x[2]))

  # Extract and name estimates
  rkhats <- c('khat' = optimization$par[1], 'rhat' = optimization$par[2])

  # End of function
  return(rkhats)
}
