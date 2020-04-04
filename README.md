# paneljudge
An R package to judge the performance of a panel of genetic markers using simulated data

## Prerequisites

The package **paneljudge** is an R package. It was developed in R version 3.5.1 (2018-07-02) using RStudio. 
To download and install R please go to [cran.r-project.org](https://cran.r-project.org).
To download and install RStudio please go to [rstudio.com](https://rstudio.com/). 
Please ensure you have the latest version of R or at least version 3.5.1. 

## Installation

A development version of **paneljudge** is available on Github. 
It can be installed in R using `install_github` from the **devtools** package.
At the time of writing (9th Oct 2019) 
an in-development version of devtools (version 2.2.1.9000)
was needed to build **paneljudge**'s vignette upon installation. 
To ensure **paneljudge** installs and the vignette builds, 
please follow the code below and accept any suggested package updates. 

```r
# Step 1) install devtools as required: 

if (!require("devtools")) { # If devtools is not intalled
  
  # Install stable version from CRAN:  
  install.packages("devtools") 
  
  # Extract and compare the stable version: 
  vdetools = as.character(packageVersion("devtools"))
  vcompare = compareVersion(vdetools, '2.2.1.9000')
  
  if (vcompare < 0) {  # If the stable version is < ‘2.2.1.9000’
    
    # Install in-development version from Github
    devtools::install_github("r-lib/devtools") 
  }

} else { # If devtools is already intalled  
  
  # Extract and compare the installed version: 
  vdetools = packageVersion("devtools")
  vcompare = compareVersion(as.character(vdetools), '2.2.1.9000')
  
  if (vcompare < 0) { # If the installed version is < ‘2.2.1.9000’
    
    # Install in-development version from Github
    devtools::install_github("r-lib/devtools") 
  }
}

# Step 2) install pixelate from GitHub 
devtools::install_github("artaylor85/pixelate", build_vignettes = TRUE, dependencies = TRUE)
```

## Usage

The intended usage of the **paneljudge** package is to judge the performance of a panel of genetic markers designed for relatedness inference. Performance is judged using data simulated under a hidden Markov model described in [1]. The package is very minimal. To see its full range of capabilities, simply load and attach **paneljudge** then read the **paneljudge** vignette accessed via XXX. Possible future additions to the package a listed below. If you would like to contribute the package see **Contributing** below or email ataylor@hsph.harvard.edu. 

## Future work

- Add Rshiny plot of marker positions with dynamic annotation: marker name, cardinality, effective cardinality, diversity etc. 
- Integrate hmmIBD [Schaffner et al. 2018, Malaria journal, 17(1)]
- Consider simulation and inference under an independent model

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)

<!--- ## Acknowledgements 
Thank you to xxxx for help testing package installation. --->
