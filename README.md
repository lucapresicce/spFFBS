# Spatiotemporal Propagation for Multivariate Bayesian Dynamic Learning

This package provides the principal parallelized implementation of forward filtering backward sampling algorithm, along with dynamic Bayesian predictive stacking to achieve exact posterior inference avoiding simulation-based approaches for multivariate spatiotemporal models.

<!--
This package provides the principal functions to perform accelerated modeling for univariate and multivariate spatial regressions. The package is used mostly within the novel working paper *"Bayesian Transfer Learning for Artificially Intelligent Geospatial Systems: A Predictive Stacking Approach" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee, 2024+)"*. To guarantee the reproducibility of scientific results, in the [Bayesian-Transfer-Learning-for-GeoAI](https://github.com/lucapresicce/Bayesian-Transfer-Learning-for-GeoAI) repository are also available all the scripts of code used for simulations, data analysis, and results presented in the Manuscript and its Supplemental material.
-->

--------------------------------------------------------------------------------
## Roadmap

| Folder | Description |
| :--- | :---: |
| `R` | contains funtions in R |
| `src` | contains function in Rcpp/C++ |

--------------------------------------------------------------------------------
## Guided installation
Since the package is not already available on CRAN (already submitted, and hopefully soon available), we use the `devtools` R package to install. Then, check for its presence on your device, otherwise install it:
```{r, echo = F, eval = F, collapse = TRUE}
if (!require(devtools)) {
  install.packages("devtools", dependencies = TRUE)
}
```
Once you have installed *devtools*, we can proceed. Let's install the `spBPS` package!
```{r}
devtools::install_github("lucapresicce/spFFBS")
```
Cool! You are ready to start, now you too could perform **_fast & feasible_** Bayesian spatiotemporal modeling!

<!--
## Tutorial for usage
-->

--------------------------------------------------------------------------------
## Contacts

| | |
| :--- | :---: |
| Author | Luca Presicce (l.presicce@campus.unimib.it) |
| Maintainer | Luca Presicce (l.presicce@campus.unimib.it) |
<!--
| Reference | [**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee (2025+) *"Bayesian Transfer Learning for Artificially Intelligent Geospatial Systems: A Predictive Stacking Approach"*  |
-->

<!--
Maintainer: l.presicce@campus.unimib.it
Reference: **Luca Presicce** and Sudipto Banerjee (2024+) *"Accelerated Meta-Kriging for massive Spatial dataset"* 
-->

