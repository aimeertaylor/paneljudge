##############################################################################
#' In this script we generate the results used in the vignette
#' paneljudge_multipanel_comparison, since it would take too long to generate
#' within the vignette itself.
##############################################################################
rm(list = ls()) # Set up
set.seed(1) # For reproducibility
devtools::load_all()

#=============================================================================
# Simulate n genotype pairs per country, r in rs, and k in ks
#=============================================================================
# Data-generating relatedness values (named s.t. simulated_genotype_pairs are
# named accordingly)
rs <- c("0.01"=0.01, "0.25"=0.25, "0.50"=0.50, "0.75"=0.75, "0.99"=0.99)

# Data-generating switch rates (named s.t. simulated_genotype_pairs are named
# accordingly)
ks <- c('1'=1,'5'=5,'10'=10,'50'=50)

# Number of pairs to per simmulate per country, r in rs and k in ks
n <- 100

# ----------------------------------------------------------------------------
# Simulate GTseq with CTSA
# ----------------------------------------------------------------------------
ds <- markers$distances # Distances
mles_CIs_GTseq <- lapply(frequencies, function(fs) {
  lapply(ks, function(k) {
    lapply(rs, function(r) {
      lapply(1:n, function(i) {
        sim_Ys <- simulate_Ys(fs, ds, k, r) # Stimulate data
        krhat <- estimate_r_and_k(fs, ds, Ys = sim_Ys) # Estimate k and r
        CIs <- compute_r_and_k_CIs(fs, ds, krhat['khat'], krhat['rhat']) # CIs
        return(cbind(krhat, CIs)) # Return mles and CIs
      })
    })
  })
})
#save(mles_CIs_GTseq, file = '../data/mles_CIs_GTseq.RData')

# ----------------------------------------------------------------------------
# Simulate GTseq without CTSA
# ----------------------------------------------------------------------------
ds <- markers[!grepl('z', markers$Amplicon_name), 'distances'] # Distances
mles_CIs_GTseq_notCTSA <- lapply(frequencies, function(fs) {
  fs <- fs[!grepl('z', markers$Amplicon_name), ]  # Extract frequencies
  lapply(ks, function(k) {
    lapply(rs, function(r) {
      lapply(1:n, function(i) {
        sim_Ys <- simulate_Ys(fs, ds, k, r) # Stimulate data
        krhat <- estimate_r_and_k(fs, ds, Ys = sim_Ys) # Estimate k and r
        CIs <- compute_r_and_k_CIs(fs, ds, krhat['khat'], krhat['rhat']) # CIs
        return(cbind(krhat, CIs)) # Return mles and CIs
      })
    })
  })
})
#save(mles_CIs_GTseq_notCTSA, file = '../data/mles_CIs_GTseq_notCTSA.RData')

# ----------------------------------------------------------------------------
# Simulate CTSA only
# ----------------------------------------------------------------------------
ds <- markers[grepl('z', markers$Amplicon_name), 'distances'] # Distances
mles_CIs_onlyCTSA <- lapply(frequencies, function(fs) {
  fs <- fs[grepl('z', markers$Amplicon_name), ] # Extract frequencies
  lapply(ks, function(k) {
    lapply(rs, function(r) {
      lapply(1:n, function(i) {
        sim_Ys <- simulate_Ys(fs, ds, k, r) # Stimulate data
        krhat <- estimate_r_and_k(fs, ds, Ys = sim_Ys) # Estimate k and r
        CIs <- compute_r_and_k_CIs(fs, ds, krhat['khat'], krhat['rhat']) # CIs
        return(cbind(krhat, CIs)) # Return mles and CIs
      })
    })
  })
})
#save(mles_CIs_onlyCTSA, file = '../data/mles_CIs_onlyCTSA.RData')

# ----------------------------------------------------------------------------
# Simulate sanger barcode
# ----------------------------------------------------------------------------
ds <- markers_sanger_barcode$distances # Distances
mles_CIs_sanger_barcode <- lapply(frequencies_sanger_barcode, function(fs) {
  lapply(ks, function(k) {
    lapply(rs, function(r) {
      lapply(1:n, function(i) {
        sim_Ys <- simulate_Ys(fs, ds, k, r) # Stimulate data
        krhat <- estimate_r_and_k(fs, ds, Ys = sim_Ys) # Estimate k and r
        CIs <- compute_r_and_k_CIs(fs, ds, krhat['khat'], krhat['rhat']) # CIs
        return(cbind(krhat, CIs)) # Return mles and CIs
      })
    })
  })
})
#save(mles_CIs_sanger_barcode, file = '../data/mles_CIs_sanger_barcode.RData')

mles_CIs_sanger_barcode$Columbia$`1`$`0.01`


