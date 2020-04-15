# paneljudge
An R package to judge the performance of a panel of genetic markers using simulated data

## Prerequisites

The package **paneljudge** is an R package. It was developed in R version 3.5.1 (2018-07-02) using RStudio. 
To download and install R please go to [cran.r-project.org](https://cran.r-project.org).
To download and install RStudio please go to [rstudio.com](https://rstudio.com/). 
Please ensure you have the latest version of R or at least version 3.5.1. 

## Installation

A development version of **paneljudge** is available on Github. 
It can be installed in R using `devtools::install_github()` from the **devtools** package (version 2.2.2).
To ensure **paneljudge** installs and the vignette builds, please follow the code below (to ensure you have **devtools** version 2.2.2 or higher) and accept any suggested package updates. 

```r
# Install or update latest stable version of devtools from CRAN
install.packages("devtools")

# Install paneljudge from GitHub 
devtools::install_github("artaylor85/paneljudge", build_vignettes = TRUE)
```

## Usage

The intended usage of the **paneljudge** package is to judge the performance of a panel of genetic markers designed for relatedness inference. Performance is judged using data (pairs of haploid genotypes) that are simulated under a hidden Markov model of relatedness [1] using allele frequency estimates provided by the user and marker positions. The markers are treated as categorical random variables whose realisations (alleles) are unordered. The effective cardinalities and diversities of the markers are computed using the input allele frequency estimates. Panel performance is judged in terms of the root mean square error (RMSE) and confidence interval width of estimated relatedness (see example figure below), where relatedness is estimated under the same model used to simulate the data. 

At present, we do not consider model misspecification; we do not account for uncertainty around input allele frequency estimates; we do not consider relatedness between pairs of haploid genotypes simulated using different allele frequencies; we do not account for marker drop-out (markers that fail to produce useful data, e.g. because they a monomorphic). Otherwise stated, at present, the performance of a panel is judged in its most favourable light; it will likely perform less well in reality. 

The package is very minimal. To see its full range of capabilities, simply load and attach **paneljudge**, consult the documentation, read the **paneljudge** vignette and view its source code (see code block below). You can also follow the example analysis stored in https://github.com/artaylor85/paneljudge/tree/master/Analysis_multipanel_multicountry/ by reading `multipanel_multicountry.pdf` and consulting the R markdown script that generated it (`multipanel_multicountry.Rmd`). 

```r
# Load and attach package
library(paneljudge)

# Lists available functions and example data sets and their documentation
help(package = "paneljudge")

# Load the vignette
vignette("paneljudge_example")

# View the all the vignette's source code, including that used to generate plots
edit(vignette("paneljudge_example"))
```

Possible future additions to the package are listed below. If you would like to contribute, see **Contributing** below or email ataylor@hsph.harvard.edu. 


![An example plot of confidence intervals around relatedness estimates based on data simulated for four different panels using frequencies from four different countries](https://github.com/artaylor85/paneljudge/blob/master/Analysis_multipanel_multicountry/multipanel_multicountry_files/figure-latex/plot%20CIs-1.pdf)


## Future work
- Impact of LD: difference between treating within microhaplotypes SNPs individually or as microhaplotypes. 
- Add tests, errors and warnings for unexpected input (frequencies etc.) 
- Add Rshiny plot of marker positions with dynamic annotations inc. marker name, effective cardinality, diversity etc. 
- Integrate hmmIBD [2].

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)

## References 
[1] Taylor et al. Genetics 212.4 (2019): 1337-1351.

[2] Schaffner et al. Malaria journal 17.1 (2018): 196.

## Acknowledgements 
Thank you to Emily LaVerriere for useful feedback. 
