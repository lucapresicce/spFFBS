# spFFBS <a href="https://lucapresicce.github.io/spFFBS/"><img src="man/figures/logo.png" alt="spFFBS website" align="right" height="138"/></a>

<!-- badges: end -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spFFBS?color=blue)](https://CRAN.R-project.org/package=spFFBS) [![Downloads](https://cranlogs.r-pkg.org/badges/spFFBS?color=orange)](https://CRAN.R-project.org/package=spFFBS) [![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/spFFBS?color=green)](https://CRAN.R-project.org/package=spFFBS)

<!-- 
# Spatiotemporal Propagation for Multivariate Bayesian Dynamic Learning
-->

## Overview

This package provides the implementation of the Forward Filtering Backward Sampling (FFBS) algorithm with Dynamic Bayesian Predictive Stacking (DYNBPS) integration for multivariate spatiotemporal models, as introduced in *"[Adaptive Markovian Spatiotemporal Transfer Learning in Multivariate Bayesian Modeling](https://arxiv.org/abs/2602.08544)" ([**Presicce**](https://lucapresicce.github.io/) and Banerjee, 2026+)*. To guarantee the reproducibility of scientific results, in the [Markovian-Spatiotemporal-Propagation](https://github.com/lucapresicce/Markovian-Spatiotemporal-Propagation) repository also includes all the scripts of code used for simulations, data analysis, and results presented in the Manuscript and its Supplemental material.

<!--
This package provides the principal functions to perform accelerated modeling for univariate and multivariate spatial regressions. The package is used mostly within the novel working paper *"Bayesian Transfer Learning for Artificially Intelligent Geospatial Systems: A Predictive Stacking Approach" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee, 2024+)"*. To guarantee the reproducibility of scientific results, in the [Bayesian-Transfer-Learning-for-GeoAI](https://github.com/lucapresicce/Bayesian-Transfer-Learning-for-GeoAI) repository are also available all the scripts of code used for simulations, data analysis, and results presented in the Manuscript and its Supplemental material.
-->

## Installation

If installing from CRAN, use the following.

``` r
install.packages("spFFBS")
```

For a quick installation of the development version, run the following command in R. We use the `devtools` R package to install. Then, check for its presence on your device, otherwise install it:

``` r
if (!require(devtools)) {
  install.packages("devtools", dependencies = TRUE)
}
```

Once you have installed *devtools*, we can proceed. Let's install the `spFFBS` package!

``` r
devtools::install_github("lucapresicce/spFFBS")
```

## Usage

Once successfully installed, load the library in R.

``` r
library(spFFBS)
```

Cool! You are ready to start, now you too could perform ***fast & feasible*** Bayesian spatiotemporal modeling!

<!--
## Tutorial for usage
-->

## Contacts

|            |                                                                                                                                                                                                                 |
|:-------------------------------|:--------------------------------------:|
| Author     |                                                           Luca Presicce ([l.presicce\@campus.unimib.it](mailto:l.presicce@campus.unimib.it))                                                            |
| Maintainer |                                                           Luca Presicce ([l.presicce\@campus.unimib.it](mailto:l.presicce@campus.unimib.it))                                                            |
| Reference  | [**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee (2026+) *"[Adaptive Markovian Spatiotemporal Transfer Learning in Multivariate Bayesian Modeling](https://arxiv.org/abs/2602.08544)"* |

<!--
Maintainer: l.presicce@campus.unimib.it
Reference: **Luca Presicce** and Sudipto Banerjee (2024+) *"Accelerated Meta-Kriging for massive Spatial dataset"* 
-->
