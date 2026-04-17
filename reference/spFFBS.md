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
  spatial = NULL,
  num_threads = 1
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

- num_threads:

  Number of cores for parallel computing (default: 1) .

## Value

A list with the components executed according to the flags.

## Examples

``` r
# \donttest{

n <- 50
t <- 5
p <- 2
q <- 2

Y <- array(rnorm(n*q*t), dim = c(n, q, t))
P <- array(rnorm(n*(p+n)*t), dim = c(n, (p+n), t))
G <- array(rnorm((p+n)*(p+n)*t), dim = c((p+n), (p+n), t))
coords <- matrix(runif(n*2), ncol = 2)
D <- as.matrix(dist(coords))

priors <- list("m" = matrix(0, n+p, q),
               "C" = diag(p+n),
               "nu" = 3,
               "Psi" = diag(q))

hyperpar <- list(tau = 0.5, phi = 1)

res <- spFFBS(Y = Y, G = G, P = P, D = D, grid = hyperpar, prior = priors)
#> 
#> ====================================================
#>          Welcome to spFFBS Bayesian Engine
#> ====================================================
#> 
#> Building parameter grid ... OK ( 1 models )
#> 
#> Running Forward Filtering (FF)...
#> 0.004 sec elapsed
#> FF completed.
#> 
#> Computing stacking weights ...
#> 
#> ====================================================
#>       Computing Dynamic Bayesian Predictive
#>              Stacking Weights (Wi)
#> ====================================================
#> 
#> Found:
#>    + Time points:   5 
#>    + Locations:     50 
#>    + Models (J):    1 
#> 
#> Using SERIAL backend (foreach sequential)
#> 
#> Evaluating model scores and optimizing weights...
#> Weight computation: 0.028 sec elapsed
#> 
#>  Weight matrix computed successfully.
#>    Dimensions:  50 x 1 
#> 
#> ====================================================
#> 
#> Total time: 0.036 sec elapsed
#> 0.037 sec elapsed
#> Global weights computed.
#> 
#> ====================================================
#>      spFFBS pipeline completed successfully!
#> ====================================================
#> 

# }
```
