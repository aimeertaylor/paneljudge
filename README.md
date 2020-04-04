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

The intended usage of the **paneljudge** package is to judge the performance of a panel of genetic markers designed for relatedness inference. Performance is judged using data simulated under a hidden Markov model described in [1]. The package is very minimal. To see its full range of capabilities, simply load and attach **paneljudge** then read the **paneljudge** vignette accessed via `vignette("paneljudge_example")`. Possible future additions to the package are listed below. If you would like to contribute, see **Contributing** below or email ataylor@hsph.harvard.edu. 

## Future work

- Add Rshiny plot of marker positions with dynamic annotations: marker name, effective cardinality, diversity etc. 
- Integrate hmmIBD [2].
- Consider simulation and inference under an independent model.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)

## References 
[1] Taylor et al. Genetics 212.4 (2019): 1337-1351.

[2] Schaffner et al. Malaria journal 17.1 (2018): 196.

<!--- ## Acknowledgements 
Thank you to xxxx for help testing package installation. --->
