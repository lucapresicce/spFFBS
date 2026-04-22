# compute_Wt: Dynamic Bayesian Predictive Stacking Weights

Computes predictive stacking weights from the forward-filter objects.
Provides a simplified and user-friendly interface with optional
parallelization.

## Usage

``` r
compute_Wt(out_FF, tau, phi, parallel = FALSE, verbose = TRUE)
```

## Arguments

- out_FF:

  Output of spFF3 (list of filtering results).

- tau:

  Vector of tau grid values (only used for column naming).

- phi:

  Vector of phi grid values (only used for column naming).

- parallel:

  Logical; use parallel backend (foreach + doParallel)? Default FALSE.

- verbose:

  Logical; print progress messages to the console? Default TRUE.

## Value

Matrix of weights of size n x J
