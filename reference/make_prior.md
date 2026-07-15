# Convert a full prior covariance matrix to block form for spFFBS

Convert a full prior covariance matrix to block form for spFFBS

## Usage

``` r
make_prior(m0, C0, nu0, Psi0, p)
```

## Arguments

- m0:

  (p+n) x q prior mean matrix.

- C0:

  (p+n) x (p+n) prior covariance matrix (full).

- nu0:

  Prior degrees of freedom (scalar).

- Psi0:

  q x q prior scale matrix.

- p:

  Number of regression coefficients (integer).

## Value

A list with elements m, C_bb, C_bw, C_ww, nu, Psi ready to pass as
`prior` to `spFFBS`.
