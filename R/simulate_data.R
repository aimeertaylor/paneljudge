#################################################################################
# function to simulate data given parameters, given fs and gendist
# fs should be a matrix with ndata rows
# and a number of columns equal to the max number of possible allele types
# if a certain position has a smaller number of possible allele types 
# then the corresponding row of fs should contain strictly positive values and then zeros.
#
# Example: if 3 types are possible at position p, but 5 types are possible
# somewhere else, then fs should have five columns, and fs[p,] might look like
# [0.2, 0.3, 0.5, 0, 0]

# gendist should be a vector where gendist[p] contains the distance between position p and p+1
# or equivalently, gendist[p-1] contains the distance between position p-1 and p, for p > 1.
# so only ndata-1 first entries of gendist are being used.
#################################################################################
simulate_data <- function(fs, gendist, k, r, epsilon = 0.001, rho = 7.4 * 10^(-7)){
  ndata <- dim(fs)[1]
  maxnstates <- dim(fs)[2]
  nstates <- 0
  Ys <- matrix(NA, nrow = ndata, ncol = 2)
  for (position in 1:ndata){
    if (position == 1){
      IBD_current <- (runif(1) <= r) # Bernoulli(r)
    } else {
      if (IBD_current){
        IBD_current <- (runif(1) < (1 - (1-r)*(1 - exp(-k * rho * gendist[position-1]))))
      } else {
        IBD_current <- (runif(1) < r*(1 - exp(-k * rho * gendist[position-1])))
      }
    }
    # number of possible allele types at current position
    nstates <- 1
    while((nstates < maxnstates) && (fs[position,nstates] > 1e-10)){
      nstates <- nstates + 1
    }
    #
    Gi <- sample(x = 0:(nstates-1), size = 1, prob = fs[position,1:nstates])
    Gj <- NA; Yi <- NA; Yj <- NA
    # generate Gi, Gj given IBD_current
    if (IBD_current){
      Gj <- Gi
    } else {
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


