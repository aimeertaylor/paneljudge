# Function to check if frequencies are ordered, sum to one, in [0,1], and informative.

# In hmmloglikelihood.cpp any frequency greater than 1e-20 is considered non-zero
fs_checks <- function(fs, non_zero_fs_lb = 1e-20, return_fs01 = FALSE){

  # Convert fs into numeric logic for zero/non-zero
  fs01 <- 1 * (fs > non_zero_fs_lb)

  # Check for disordered frequencies.
  # Specifically, check fs ordered s.t. non-zero precede zero for all markers.
  # Disordered frequencies break the hmmloglikelihood.cpp code.
  if(!all(apply(fs01, 1, function(x) all(x == cummin(x))))) {
    stop ("Disordered fs. Per row, all non-zero frequencies should precede all zero frequencies.")
  }

  # Check per-marker frequencies sum to one
  # Those that do not do not break the hmmloglikelihood.cpp code so are liable to go undetected
  if(any(rowSums(fs) < 1-non_zero_fs_lb) | any(rowSums(fs) > 1 + non_zero_fs_lb)) {
    stop("Some markers have frequencies that don't sum to one.")
  }

  # Check frequencies are in [0,1]
  # Those that exceed one do not break the hmmloglikelihood.cpp code so are liable to go undetected
  # Those that are negative break the hmmloglikelihood.cpp code
  if(any(fs < 0) | any(fs > 1)) {
    stop("Some frequencies are not in [0,1].")
  }

  # Check for uninformative markers
  # These markers do not break the hmmloglikelihood.cpp code but could be omitted.
  if (any(fs == 1)) {
    warning("Some markers are uninformative (have allele frequencies equal to one).")
  }

  if(return_fs01) return(fs01) # For use in Ys check
}






