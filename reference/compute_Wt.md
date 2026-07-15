# Compute dynamic Bayesian predictive stacking weights (v2)

Calls `compute_Wt_cpp_v2()` directly - no foreach overhead.

## Usage

``` r
compute_Wt(out_FF, tau, phi, n_threads = 1L, verbose = TRUE)
```

## Arguments

- out_FF:

  Output of `spFF3_v2`.

- tau:

  Numeric vector of tau grid values (for column names only).

- phi:

  Numeric vector of phi grid values (for column names only).

- n_threads:

  Number of OpenMP threads (default 1).

- verbose:

  Logical; print progress? Default TRUE.

## Value

n x J weight matrix.
