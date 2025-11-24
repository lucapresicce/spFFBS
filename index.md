# Spatiotemporal Propagation for Multivariate Bayesian Dynamic Learning

This package provides the principal parallelized implementation of
forward filtering backward sampling algorithm, along with dynamic
Bayesian predictive stacking to achieve exact posterior inference
avoiding simulation-based approaches for multivariate spatiotemporal
models.

| \## Roadmap                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| \## Guided installation Since the package is not already available on CRAN (already submitted, and hopefully soon available), we use the `devtools` R package to install. Then, check for its presence on your device, otherwise install it: `{r, echo = F, eval = F, collapse = TRUE} if (!require(devtools)) { install.packages("devtools", dependencies = TRUE) }` Once you have installed *devtools*, we can proceed. Let’s install the `spBPS` package! `{r} devtools::install_github("lucapresicce/spFFBS")` Cool! You are ready to start, now you too could perform ***fast & feasible*** Bayesian spatiotemporal modeling! |
|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |

## Contacts

|            |                                               |
|:-----------|:---------------------------------------------:|
| Author     | Luca Presicce (<l.presicce@campus.unimib.it>) |
| Maintainer | Luca Presicce (<l.presicce@campus.unimib.it>) |
