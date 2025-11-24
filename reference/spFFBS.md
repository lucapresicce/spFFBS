# spFFBS: Spatiotemporal Bayesian Pipeline (friendly interface)

A user-friendly, modular wrapper for running a Bayesian spatiotemporal
filtering + weighting pipeline, with optional backward sampling,
forecasting, and spatial interpolation.

## Usage

``` r
spFFBS(
  Y,
  G,
  P,
  D,
  grid = list(tau = NULL, phi = NULL),
  prior,
  do_BS = FALSE,
  do_forecast = FALSE,
  do_spatial = FALSE,
  L = 200,
  tnew = NULL,
  spatial = NULL
)
```

## Arguments

- Y:

  Response data (3D array or cube).

- G:

  System matrix (cube).

- P:

  Observation matrix (cube).

- D:

  Spatial distance matrix.

- grid:

  List with elements:

  - tau: numeric vector

  - phi: numeric vector

- prior:

  Prior list for forward filter (m, C, nu, Psi)

- do_BS:

  Logical: run backward sampling? (default: FALSE)

- do_forecast:

  Logical: run temporal forecasts? (default: FALSE)

- do_spatial:

  Logical: run spatial interpolation? (default: FALSE)

- L:

  Number of posterior samples (default 200)

- tnew:

  Forecast horizon (default 5)

- spatial:

  Optional list for spatial: list(crd = , crdtilde = , Xtilde = )

## Value

A list with the components executed according to the flags.
