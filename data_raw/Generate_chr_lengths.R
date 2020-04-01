# Chromosome lengths (in base pairs) copied from plasmoDB (there are better ways)
chr_lengths <- c("1"=640851,
                "2"=947102,
                "3"=1067971,
                "4"=1200490,
                "5"=1343557,
                "6"=1418242,
                "7"=1445207,
                "8"=1472805,
                "9"=1541735,
                "10"=1687656,
                "11"=2038340,
                "12"=2271494,
                "13"=2925236,
                "14"=3291936)

save(chr_lengths, file = '../data/chr_lengths.RData')
