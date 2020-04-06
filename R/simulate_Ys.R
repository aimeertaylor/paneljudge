###########################################################################
#' Simulate genotype calls for a pair of haploid genotypes
#'
#' Given a matrix of marker allele frequencies, a vector of inter-marker
#' distances, a relatedness parameter, and a switch rate parameter, for a pair
#' of haploid genotypes \code{simulate_Ys} returns genotype calls simulated
#' under the HMM described in [1].
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
#' @param k data-generating switch rate parameter, i.e. \eqn{k} in [1].
#' @param r data-generating relatedness parameter, i.e. \eqn{r} in [1].
#' @param epsilon Genotyping error, i.e. \eqn{\epsilon} in [1]. The genotyping
#'   error is the probability of miscalling one specific allele for another. As
#'   such, the error rate for the t-th marker, \eqn{(Kt-1)\epsilon}, scales with
#'   \eqn{Kt} (the per-marker allele count, cardinality).
#' @param rho Recombination rate, i.e. \eqn{\rho} in [1]. The recombination rate
#'   corresponds to the probability of a crossover per base pair. It is assumed
#'   constant across the genome under the HMM of [1]. Its default value
#'   corresponds to an average rate estimated for \emph{Plasmodium falciparum}
#'   [2].
#'
#' @return Simulated genotype calls for a pair of haploid genotypes, i.e. the
#'   \eqn{Yt}s of the \eqn{i}-th and \eqn{j}-th haploid genotypes in [1].
#'   Specifically, a \eqn{m} by 2 matrix, where \eqn{m} is the marker count and
#'   each column contains a haploid genotype. For all \eqn{t = 1,...,m} markers,
#'   alleles are enumerated 0 to \eqn{Kt-1}, where \eqn{Kt} is the cardinality
#'   (per-marker allele count) of the \eqn{t}-th marker. For example, if \eqn{Kt
#'   = 2}, both \code{Ys[t,1]} and \code{Ys[t,2]} are either 0 or 1.
#'
#' @examples
#' simulate_Ys(fs = frequencies$Colombia, ds = markers$distances, k = 10, r = 0.5)
#'
#' @references \enumerate{ \item Taylor, A.R., Jacob, P.E., Neafsey, D.E. and
#'   Buckee, C.O., 2019. Estimating relatedness between malaria parasites.
#'   Genetics, 212(4), pp.1337-1351. \item Miles, A., Iqbal, Z., Vauterin, P.,
#'   Pearson, R., Campino, S., Theron, M., Gould, K., Mead, D., Drury, E.,
#'   O'Brien, J. and Rubio, V.R., 2016. Indels, structural variation, and
#'   recombination drive genomic diversity in Plasmodium falciparum. Genome
#'   research, 26(9), pp.1288-1299.}
#' @export
###########################################################################

simulate_Ys <- function(fs, ds, k, r, epsilon = 0.001, rho = 7.4 * 10^(-7)){

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


