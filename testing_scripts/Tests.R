# To convert into unit tests

# # Remove old
remove.packages("paneljudge")

# # Install new
devtools::install_github("artaylor85/paneljudge", build_vignettes = TRUE)

# Attach
library(paneljudge)

# Help
help(package = "paneljudge")

# Function to test
test_all <- function(fs, ds, fs_warn = TRUE, k = 5, r = 0.25, epsilon_ = 0.01) {

  # Compute diversities and effective cardinalities
  divs <- compute_diversities(fs, fs_warn)
  eff_card <- compute_eff_cardinalities(fs, fs_warn)

  # Stimulate some data
  simulated_Ys <- simulate_Ys(fs, ds, k, r, warn_fs = fs_warn, epsilon = epsilon_)

  # Estimate the switch rate parameter, k, and relatedness parameter, r
  krhat <- estimate_r_and_k(fs, ds, Ys = simulated_Ys, warn_fs = fs_warn, epsilon = epsilon_)

  # Compute confidence intervals (CIs)
  compute_r_and_k_CIs(fs, ds, khat = krhat['khat'], rhat = krhat['rhat'], warn_fs = fs_warn, epsilon = epsilon_)
}


fs <- frequencies$Colombia
ds <- markers$distances

# Default:
set.seed(1)
test_all(fs, ds, k = 5, r = 0.25)

# Warnings off:
set.seed(1)
test_all(fs, ds, k = 5, r = 0.25, fs_warn = FALSE)

# Erroneous fs: disordered rows
fs <- frequencies$Colombia
head(fs[1,]); fs[1,] <- fs[1,]+0.3; head(fs[1,])
set.seed(1)
test_all(fs, ds, k = 5, r = 0.25, fs_warn = FALSE)

# Erroneous fs: sum to one deviates by a large amount
fs <- frequencies$Colombia
fs[1,] <- fs[1,]+0.3; head(fs[,1:10])
set.seed(1)
test_all(fs, ds, k = 5, r = 0.25, fs_warn = FALSE)

# Erroneous fs: sum to one deviates by a small amount with a frequency exceeding one
fs <- frequencies$Colombia
fs[2,1] <- 1+1e-6; head(fs[,1:10])
set.seed(1)
test_all(fs, ds, k = 5, r = 0.25, fs_warn = TRUE)
test_all(fs, ds, k = 5, r = 0.25, fs_warn = FALSE)

# Erroneous fs: negative frequency
fs <- frequencies$Colombia
fs[1,1] <- -fs[1,1]; head(fs[,1:3])
set.seed(1)
test_all(fs, ds, k = 5, r = 0.25, fs_warn = FALSE)

# Epsilon check
fs <- frequencies$Colombia
set.seed(1)
# Stimulate some data with lots of error then try to evaluate with epsilon = 0
simulated_Ys <- simulate_Ys(fs, ds, k = 5, r = 0.25, warn_fs = FALSE, epsilon = 0.5)
krhat <- estimate_r_and_k(fs, ds, Ys = simulated_Ys, warn_fs = FALSE, epsilon = 0)

# NAs in Ys check
simulated_Ys[1,1] <- NA
krhat <- estimate_r_and_k(fs, ds, Ys = simulated_Ys, warn_fs = FALSE, epsilon = 0)


set.seed(1)
fs <- frequencies$Colombia
test_all(fs, ds, k = -5, r = 0.25, fs_warn = FALSE)
test_all(fs, ds, k = 5, r = -0.25, fs_warn = FALSE)
test_all(fs, ds, k = -5, r = 0.25, fs_warn = FALSE)
test_all(fs, ds, k = 5, r = 1.25, fs_warn = FALSE)


