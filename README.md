# drought4R

A drought assessment package for the [climate4R](http://www.meteo.unican.es/climate4r) Framework for climate data access and analysis.

This package is a wrapper to the R package SPEI (Beguería and Serrano 2017), tailored to its usage with large climate model datasets within climate4R.

## Installation

The latest stable release can be installed using `install_github` from package `devtools` (Wickham _et al._ 2018):

```R
devtools::install_github("SantanderMetGroup/drought4R")
```
## Worked examples

A worked example showing how to compute bias-corrected SPEI projections using Euro-CORDEX data from the User Data Gateway is provided in the package documentation. Once the package has been installed, type:

```R
utils::RShowDoc("drought4R_notebook", package = "drought4R")
```

    

