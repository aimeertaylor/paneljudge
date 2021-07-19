# Function to check if frequencies are ordered, sum to one, in [0,1], and informative.

# In hmmloglikelihood.cpp any frequency greater than 1e-20 is considered non-zero
fs_checks <- function(fs, warn,
                      non_zero_fs_lb = 1e-20, # Magic number
                      max_dev_lim = 1e-5, # Magic number
                      do_return = FALSE){

  if (!warn) {
    # turn off new warn option while returning old
    op <- options(warn = -1)
    on.exit(options(op))
  }

  # Convert fs into numeric logic for zero/non-zero
  fs01 <- 1 * (fs > non_zero_fs_lb)

  # Check frequencies are in [0,1]
  # Those that exceed one do not break the hmmloglikelihood.cpp code so are liable to go undetected
  # Those that are negative break the hmmloglikelihood.cpp code
  # Those that are negative also trigger error below due to cummin so this check precedes one below
  if(any(fs < 0) | any(fs > 1)) {
    stop("Some frequencies are not in [0,1].", call. = FALSE)
  }

  # Check for disordered frequencies.
  # Specifically, check fs ordered s.t. non-zero precede zero for all markers.
  # Disordered frequencies break the hmmloglikelihood.cpp code.
  if(!all(apply(fs01, 1, function(x) all(x == cummin(x))))) {
    stop ("Disordered fs. Per row, all non-zero frequencies should precede all zero frequencies.",
          call. = FALSE)
  }

  # Check per-marker frequencies sum to one
  # Those that do not do not break the hmmloglikelihood.cpp code so are liable to go undetected
  fs_sum <- rowSums(fs)
  if(any(fs_sum != 1)) {
    max_dev <- max(abs(fs_sum - 1))
    max_dev
    if(max_dev > max_dev_lim) {
      stop(sprintf("Some markers have frequencies whose sum deviates from one by up to %s.", max_dev),
           call. = FALSE)
    } else {
      warning(sprintf("Some markers have frequencies whose sum deviates from one by up to %s.", max_dev),
              call. = FALSE)
    }
  }



  # Check for uninformative markers
  # These markers do not break the hmmloglikelihood.cpp code but could be omitted.
  if (any(fs == 1)) {
    warning("Some markers are uninformative (have allele frequencies equal to one).",
            call. = FALSE)
  }

  if(do_return) { # For use in Ys check in estimate_r_and_k()
    return(list(fs01 = fs01,
                non_zero_fs_lb = non_zero_fs_lb))
    }
}






