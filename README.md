<!-- README.md is generated from README.Rmd. Please edit that file -->




# mi4p <img src="man/figures/logo.png" align="right" width="200"/>

# A multiple imputation framework for proteomics
## Marie Chion, Christine Carapito and Frédéric Bertrand

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/mariechion/mi4p/workflows/R-CMD-check/badge.svg)](https://github.com/mariechion/mi4p/actions)
<!--[![Codecov test coverage](https://codecov.io/gh/mariechion/mi4p/branch/master/graph/badge.svg)](https://codecov.io/gh/mariechion/mi4p?branch=master)-->
[![CRAN status](https://www.r-pkg.org/badges/version/mi4p)](https://cran.r-project.org/package=mi4p)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/mi4p)](https://cran.r-project.org/package=mi4p)
[![GitHub Repo stars](https://img.shields.io/github/stars/mariechion/mi4p?style=social)](https://github.com/mariechion/mi4p)
<!--[![DOI](https://zenodo.org/badge/18441799.svg)](https://zenodo.org/badge/latestdoi/18441799)-->
<!-- badges: end -->

This repository contains the R code and package for the _mi4p_ methodology (Multiple Imputation for Proteomics), proposed by Marie Chion, Christine Carapito and Frédéric Bertrand. 

The following material is available on thr Gthub repository of the package (https://github.com/mariechion/mi4p/).

1. The `Functions` folder contains all the functions used for the workflow.  

2. The `Simulation-1`, `Simulation-2` and `Simulation-3` folders contain all the R scripts and data used to conduct simulated experiments and evaluate our methodology. 

3. The  `Arabidopsis_UPS` and `Yeast_UPS` folders contain all the R scripts and data used to challenge our methodology on real proteomics datasets. Raw experimental data were deposited with the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifiers PXD003841 and PXD027800.

This website and these examples were created by M. Chion, C. Carapito and F. Bertrand.

## Installation

You can install the released version of mi4p from [CRAN](https://CRAN.R-project.org) with:


```r
install.packages("mi4p")
```

You can install the development version of mi4p from [github](https://github.com) with:


```r
devtools::install_github("mariechion/mi4p")
```

## Examples

### First section


```r
library(mi4p)
```

