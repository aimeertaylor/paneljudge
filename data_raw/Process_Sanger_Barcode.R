###################################################################
# This script processes data on the Sanger amplicon barcode for 
# Senegal, Columbia, and French Guiana 
# as processed by Tim Straub June 4th 2019
#
# File path and files 
#/gsap/garage-protistvector/tstraub/ampseq/Sanger/for_Aimee
#
# All_Sanger_Amplicon_Haplotype_Frqs.Rdata          - haplotype frequencies for all three
# Columbia_Sanger_Amplicon_Haplotype_Frqs.Rdata	    - haplotype frequencies for Colombia
# FrenchGuiana_Sanger_Amplicon_Haplotype_Frqs.Rdata - haplotype frequencies for French Guiana
# Senegal_Sanger_Amplicon_Haplotype_Frqs.Rdata      - haplotype frequencies for Senegal	  
# combined_amp_frqs.txt                             - haplotype frequencies?
# combined_amp_frqs_sparse_matrix.txt               - haplotype frequencies?
#
# scp aimeet@login.broadinstitute.org:/gsap/garage-protistvector/tstraub/ampseq/Sanger/for_Aimee/* ./RData/
###################################################################
rm(list = ls())

# Load data
load('../RData/All_Sanger_Amplicon_Haplotype_Frqs.RData')

# Convert factors to character strings 
amp_frqs[] <- lapply(amp_frqs, as.character)

# First separate amplicon data from country specific data 
amp_data = amp_frqs[!duplicated(amp_frqs$Amplicon_name), c('Amplicon_name','Chr','Start','Stop')]
rownames(amp_data) = NULL

# Add median positions
amp_data$pos = apply(amp_data[,c('Start','Stop')], 1, function(x)median(as.numeric(x)))

# Add numeric chromosome
amp_data$chrom = as.numeric(gsub('_v3', '', gsub('Pf3D7_', '', amp_data$Chr)))
table(amp_data$chrom) # Check all chromosomes are represented: yes 

# Re-order data 
reordered_data_list = list()
for(chr in sort(unique(amp_data$chrom))){
  x = amp_data$pos[chr == amp_data$chrom] # Extract positions per chromosome 
  print(all(x == cummax(x))) # If not all True, not all are monotonically increasing
  inds = sort.int(x, index.return = T)$ix # Sort positions
  reordered_data_list[[chr]] = amp_data[chr == amp_data$chrom, ][inds, ]
}

# Re-ordered data
amp_data = do.call(rbind, reordered_data_list)
all(sapply(unique(amp_data$chrom), function(chr){ # Check order
  x = amp_data$pos[chr == amp_data$chrom] 
  all(x == cummax(x)) 
}))

# Create distances
amp_data$dt <- c(diff(amp_data$pos), Inf)
pos_change_chrom <- 1 + which(diff(amp_data$chrom) != 0) # find places where chromosome changes
amp_data$dt[pos_change_chrom-1] <- Inf



# ************************ IMPORTANT ***************************
# After this next step, due to ordering of frequencies, the allele of
# "Genotype.1" in Amplicon #1 in Sengal may no longer be the 
# same allele as that of "Genotype.1" in Amplicon #1 in Colombia
# **************************************************************

# Separate by Country, order amplicons, order frequencies 
Genotypes = names(amp_frqs)[grepl('Genotype', names(amp_frqs))]
amp_order = amp_data$Amplicon_name
frqs_per_country = dlply(amp_frqs, 'Country', function(x){
  rownames(x) = x$Amplicon_name
  y = x[amp_order,] # Order amplicons
  z = t(apply(y[, Genotypes], 1, function(f) sort(as.numeric(f), decreasing = T))) # Order frequencies as numeric
})


# Save for reference
save(amp_data,frqs_per_country, file = '../RData/sanger_amp_data_for_IBDsim.RData')



load('../RData/sanger_amp_data_for_IBDsim.RData')

Cardinality_summaries = sapply(names(frqs_per_country), function(Country){
  Kts = 1/rowSums(frqs_per_country[[Country]]^2)
  range_Kt = range(Kts)
  c(mean = mean(Kts), 
    min = range_Kt[1], 
    max = range_Kt[2], 
    stdev = sd(Kts))
})

round(Cardinality_summaries, 2)
