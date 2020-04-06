###########################################################################
#' Simulate data on a pair of haploid genotypes
#'
#' For a vector of input frequencies and distances, \code{simulate_data()}
#' returns simulated data on a pair of haploid genotypes.
#'
#' @param fs frequencies. m (marker count) by Kmax (max cardinality over m
#'   markers) matrix if Kt < Kmax for any t in 1:m, then fs[t,1:Kt] in (0,1) &
#'   fs[t,(Kt+1):Kmax] = 0 Example: if Kt = 2 < Kmax = 4 then fs[t,] might look
#'   like [0.2, 0.7, 0, 0].
#' @param ds distances vector where ds[t] contains the distance
#'   between t t and t+1 or equivalently, ds[t-1] contains the
#'   distance between t t-1 and t, for p > 1. s.t. only m-1 first entries
#'   of ds are being used.
#' @param k data-generating switch rate parameter.
#' @param r data-generating relatedness parameter.
#' @param epsilon genotyping error (probability of miscall).
#' @param rho recombination rate in probability of break point per base pair.
#' @return Ys simulated dated.
#' @examples
#' simulate_data(fs = frequencies$Colombia, ds = markers$dt, k = 10, r = 0.5)
#'
#' @references Taylor, A.R., Jacob, P.E., Neafsey, D.E. and Buckee, C.O., 2019.
#'   Estimating relatedness between malaria parasites. Genetics, 212(4),
#'   pp.1337-1351.
#' @export
###########################################################################

simulate_data <- function(fs, ds, k, r, epsilon = 0.001, rho = 7.4 * 10^(-7)){

  m <- dim(fs)[1] # Extract marker count
  Kmax <- dim(fs)[2] # Extract Kmax
  Ys <- matrix(NA, nrow = m, ncol = 2) # Create simulated data store

  for(t in 1:m){ # For t = 1,...,m

    if (t == 1){ # At t = 1, draw IBDt from Bernoulli(r)
      IBD_t <- (runif(1) <= r)
    } else { # For t > 1, draw IBDt according to the transition matrix
      if (IBD_t){
        IBD_t <- (runif(1) < (1 - (1-r)*(1 - exp(-k * rho * ds[t-1]))))
      } else {
        IBD_t <- (runif(1) < r*(1 - exp(-k * rho * ds[t-1])))
      }
    }

    # Extract the t-th marker cardinality (allele count)
    Kt <- 1
    while((Kt < Kmax) && (fs[t,Kt] > 1e-10)){
      Kt <- Kt + 1
    }

    # Sample Gi (true alleles of the i-th individual) according to the allele frequencies
    Gi <- sample(x = 0:(Kt-1), size = 1, prob = fs[t,1:Kt])

    # Sample Gj given IBD_t and Gi
    if (IBD_t){ # Copy true alleles of Gi
      Gj <- Gi
    } else { # Draw independently according to the allele frequencies
      Gj <- sample(x = 0:(Kt-1), size = 1, prob = fs[t,1:Kt])
    }

    # Generate Yi given Gi (i.e. potentially erroneous genotyping call)
    if (runif(1) < (1 - (Kt-1)*epsilon)){ # No error
      Yi <- Gi
    } else {  # Sample uniformly on the other possible alleles
      otheralleles <- setdiff(0:(Kt-1), Gi)
      Yi <- sample(x = otheralleles, size = 1)
    }

    # Generate Yj given Gj (i.e. potentially erroneous genotyping call)
    if (runif(1) < (1 - (Kt-1)*epsilon)){ # No error
      Yj <- Gj
    } else { # Sample uniformly on the other possible alleles
      otheralleles <- setdiff(0:(Kt-1), Gj)
      Yj <- sample(x = otheralleles, size = 1)
    }

    # Store simulated data
    Ys[t,] <- c(Yi,Yj)
  }

  # End of function
  return(Ys)
}


