#################################################################################
# Function to simulate data under the HMM that features in
# Taylor, A.R., Jacob, P.E., Neafsey, D.E. and Buckee, C.O., 2019.
# Estimating relatedness between malaria parasites.
# Genetics, 212(4), pp.1337-1351.
#
# fs: ndata (marker count) by Kmax (max cardinality over m markers) matrix
# if Kt < Kmax for any t in 1:m, then fs[t,1:Kt] in (0,1) & fs[t,(Kt+1):Kmax] = 0
# Example: if Kt = 2 < Kmax = 4 then fs[t,] might look like [0.2, 0.7, 0, 0]
#
# gendist: vector where gendist[t] contains the distance between position t and t+1
# or equivalently, gendist[t-1] contains the distance between position t-1 and t, for p > 1.
# s.t. only m-1 first entries of gendist are being used.
#################################################################################
simulate_data <- function(fs, gendist, k, r, epsilon = 0.001, rho = 7.4 * 10^(-7)){

  ndata <- dim(fs)[1] # Extract marker count
  maxnstates <- dim(fs)[2] # Extract Kmax
  nstates <- 0
  Ys <- matrix(NA, nrow = ndata, ncol = 2) # Create store for simulated data
  for (position in 1:ndata){ # For t = 1,...,m
    if (position == 1){ # At t = 1, draw IBDt from Bernoulli(r)
      IBD_current <- (runif(1) <= r) # draw Bernoulli(r)
    } else { # For t > 1, draw IBDt according to the transition matrix
      if (IBD_current){
        IBD_current <- (runif(1) < (1 - (1-r)*(1 - exp(-k * rho * gendist[position-1]))))
      } else {
        IBD_current <- (runif(1) < r*(1 - exp(-k * rho * gendist[position-1])))
      }
    }
    # Extract Kt number of possible allele types at current position
    nstates <- 1
    while((nstates < maxnstates) && (fs[position,nstates] > 1e-10)){
      nstates <- nstates + 1
    }
    # Sample true alleles of the i-th individual according to their frequencies
    Gi <- sample(x = 0:(nstates-1), size = 1, prob = fs[position,1:nstates])
    Gj <- NA; Yi <- NA; Yj <- NA
    # generate Gi, Gj given IBD_current
    if (IBD_current){ # Copy true allels of Gi
      Gj <- Gi
    } else { # Independtly draw
      Gj <- sample(x = 0:(nstates-1), size = 1, prob = fs[position,1:nstates])
    }
    # generate Yi, Yj given Gi, Gj,
    # i.e. genotyping error part of the model
    if (runif(1) < (1 - (nstates-1)*epsilon)){
      Yi <- Gi # no error
    } else {
      # sample uniformly on the other possible states
      otherstates <- setdiff(0:(nstates-1), Gi)
      Yi <- sample(x = otherstates, size = 1)
    }
    if (runif(1) < (1 - (nstates-1)*epsilon)){
      Yj <- Gj
    } else {
      otherstates <- setdiff(0:(nstates-1), Gj)
      Yj <- sample(x = otherstates, size = 1)
    }
    Ys[position,] <- c(Yi,Yj)
  }
  return(Ys)
}


