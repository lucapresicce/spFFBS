## Submission spFFBS (version 0.0-1)

### Package description

Provides functions for the Forward Filtering Backward Sampling (FFBS) algorithm with Dynamic Bayesian Predictive Stacking (DYNBPS) integration for multivariate spatiotemporal models, as introduced in Presicce and Banerjee (2026+) <doi:10.48550/arXiv.2602.08544>. The core functions leverage 'C++' and OpenMP for high-performance parallel computation.

### Dependencies

Dependencies kept to a minimum for portability. The package introduces: 
  - `spFFBS()`: single unified function incorporating posterior filtering, forecasting, and spatial predictive sampling to
  simplify workflow execution.

  - Working examples and a vignette are provided.

## R CMD check results

0 errors \| 0 warnings \| 1 note

**NOTE:** (checking for future file timestamps ... NOTE / unable to verify current time)

### Test environments

-   R-hub: linux (ubuntu-latest), macos-arm64 (macos-latest), windows (windows-latest)
-   Local: Windows 11, R 4.5.3

**NOTE (Win-Builder only):** `check_win_devel()` fails with a compilation error in Rcpp's own headers (`Rcpp/Function.h:75: R_NamespaceRegistry not declared`), probably caused by an outdated Rcpp in Win-Builder's `D:/RCompile/CRANpkg/lib/4.6/`. This error is not in the package code. The package compiles and checks cleanly on: 
  - R-hub: linux, macos-arm64, windows (R-release) 
  - Local Windows (R 4.5.3)
