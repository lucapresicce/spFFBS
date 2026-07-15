utils::globalVariables("i")

# ============================================================
#  FFBS-DYNBPS-v2.R
#  Optimised R wrapper for the spFFBS C++ backend.
#
#  Key differences from FFBS-DYNBPS.R (v1):
#
#  [I-1]  G accepted as a pxp matrix (G_beta) + scalar rho
#         instead of a full (p+n)x(p+n)xT cube.
#         Default G_beta = I_p, rho = 1 = identity system matrix.
#         Users wanting a custom block-diagonal G pass G_beta and rho.
#
#  [I-2]  P replaced by X: an nxpxT array of covariates.
#         The structure P_t = [X_t | I_n] is exploited algebraically
#         inside C++; the nx(p+n) matrix is never formed.
#         Eliminates up to 18 GB for n=1500, T=1000.
#
#  [I-3]  Prior now stores covariance in block form:
#         list(m, C_bb (pxp), C_bw (pxn), C_ww (nxn), nu, Psi).
#         Helper make_prior() builds this from a full C0 matrix
#         for backward compatibility.
#
#  [I-4]  compute_Wt() calls compute_Wt_cpp_v2() directly,
#         bypassing the foreach/doParallel overhead of v1.
#
#  [I-5]  Seed propagation: call seed_rng_v2() before any C++ work
#         so that set.seed() in R is respected inside parallel loops.
#
#  [I-6]  spFFBS() is fully backward-compatible in outputs:
#         the returned list has the same names as spFFBS() v1
#         (FF, Wi, Wglobal, BS, forecast, spatial).
# ============================================================


# -- Internal: propagate R seed into the C++ RNG pool ------------------------
.seed_cpp_rng <- function() {
  seed_rng_v2(stats::runif(1), stats::runif(1))
}


# -- Helper: build block prior from a full C0 matrix -------------------------

#' Convert a full prior covariance matrix to block form for spFFBS
#'
#' @param m0   (p+n) x q prior mean matrix.
#' @param C0   (p+n) x (p+n) prior covariance matrix (full).
#' @param nu0  Prior degrees of freedom (scalar).
#' @param Psi0 q x q prior scale matrix.
#' @param p    Number of regression coefficients (integer).
#'
#' @return A list with elements m, C_bb, C_bw, C_ww, nu, Psi
#'   ready to pass as \code{prior} to \code{spFFBS}.
#'
#' @export
make_prior <- function(m0, C0, nu0, Psi0, p) {
  n <- nrow(C0) - p
  if (n <= 0) stop("p must be less than nrow(C0).")
  list(
    m    = m0,
    C_bb = C0[seq_len(p),          seq_len(p),          drop = FALSE],
    C_bw = C0[seq_len(p),          p + seq_len(n),      drop = FALSE],
    C_ww = C0[p + seq_len(n),      p + seq_len(n),      drop = FALSE],
    nu   = nu0,
    Psi  = Psi0
  )
}


# -- Internal: compute stacking weights via C++ -------------------------------

#' Compute dynamic Bayesian predictive stacking weights (v2)
#'
#' Calls \code{compute_Wt_cpp_v2()} directly - no foreach overhead.
#'
#' @param out_FF  Output of \code{spFF3_v2}.
#' @param tau     Numeric vector of tau grid values (for column names only).
#' @param phi     Numeric vector of phi grid values (for column names only).
#' @param n_threads Number of OpenMP threads (default 1).
#' @param verbose Logical; print progress? Default TRUE.
#'
#' @return n x J weight matrix.
#' @keywords internal
compute_Wt <- function(out_FF, tau, phi,
                          n_threads = 1L, verbose = TRUE) {

  if (verbose) {
    cat("\n====================================================\n")
    cat("      Computing Dynamic Bayesian Predictive\n")
    cat("             Stacking Weights (Wi) - v2\n")
    cat("====================================================\n\n")
  }

  tmax <- length(out_FF)
  if (tmax < 2L) stop("Need at least two time slices for predictive stacking.")

  density0 <- out_FF[[2L]]$density_evaluations
  if (is.null(density0))
    stop("Forward filter output does not contain $density_evaluations.")

  n <- nrow(density0)
  J <- ncol(density0)

  if (verbose) {
    cat(sprintf("  Time points : %d\n  Locations   : %d\n  Models (J)  : %d\n\n",
                tmax, n, J))
    cat("  Computing weights via compute_Wt_cpp_v2 ...\n")
  }

  # Build density list (log-scale) - indices 1..(tmax-1) for t=2..tmax
  density_list <- lapply(seq(2L, tmax), function(s)
    out_FF[[s]]$density_evaluations)

  t0 <- proc.time()["elapsed"]

  # [I-4] Direct C++ call - no R-level loop, no foreach
  Wi <- compute_Wt_cpp_v2(
    density_list = density_list,
    n       = n,
    t       = tmax,
    lr      = 0.05,
    max_iter = 500L,
    n_threads = as.integer(n_threads)
  )

  elapsed <- proc.time()["elapsed"] - t0

  # Name columns
  par_grid   <- spBPS:::expand_grid_cpp(rev(tau), rev(phi))
  colnames(Wi) <- paste0("M", seq_len(nrow(par_grid)))
  attr(Wi, "model_labels") <- sprintf(
    "Model %d [ tau = %.4f, phi = %.4f ]",
    seq_len(nrow(par_grid)), par_grid[, 1L], par_grid[, 2L])

  if (verbose) {
    cat(sprintf("  Weight matrix: %d x %d  (%.2f s)\n\n", nrow(Wi), ncol(Wi), elapsed))
    cat("====================================================\n\n")
  }

  Wi
}


# -- Main user-facing function ------------------------------------------------

#' spFFBS: Spatiotemporal Bayesian Pipeline (optimised)
#'
#' @description
#' Optimised wrapper for the spFFBS v2 C++ backend. Runs forward
#' filtering, dynamic stacking weight computation, and optionally backward
#' sampling, temporal forecasting, and spatial interpolation.
#'
#' @section Differences from \code{spFFBS} (v1):
#' \itemize{
#'   \item \code{G} is now a \emph{pxp matrix} (\code{G_beta}) plus a scalar
#'     \code{rho}, representing the block-diagonal system matrix
#'     \eqn{G = \mathrm{diag}(G_{\beta},\, \rho I_n)}.
#'     Defaults to the identity (\code{G_beta = NULL, rho = 1}).
#'   \item \code{X} (nxpxT array of covariates) replaces the \code{P} cube.
#'     The observation matrix \eqn{P_t = [X_t \mid I_n]} is formed implicitly
#'     inside C++.
#'   \item \code{prior} should be in block form; use \code{make_prior()}
#'     to convert a full \eqn{C_0} matrix if needed.
#' }
#'
#' @param Y     Response array, n x q x T.
#' @param X     Covariate array, n x p x T. Replaces \code{P}.
#' @param D     n x n spatial distance matrix.
#' @param grid  List with elements \code{tau} and \code{phi} (numeric vectors).
#' @param prior Block-form prior list; see \code{\link{make_prior}}.
#'   Elements: \code{m} ((p+n)xq), \code{C_bb} (pxp), \code{C_bw} (pxn),
#'   \code{C_ww} (nxn), \code{nu} (scalar), \code{Psi} (qxq).
#'   If the old-style full \code{C} is supplied instead of \code{C_bb/C_bw/C_ww},
#'   the wrapper auto-converts using \code{make_prior()}.
#' @param G_beta  pxp upper-left block of the system matrix G.
#'   Default \code{NULL} = identity.
#' @param rho     Scalar multiplier for the lower-right nxn identity block of G.
#'   Default 1.
#' @param do_BS      Logical; run backward sampling? Default FALSE.
#' @param do_forecast Logical; run temporal forecast? Default FALSE.
#' @param do_spatial  Logical; run spatial interpolation? Default FALSE.
#' @param L        Number of posterior samples. Default 200.
#' @param tnew     Forecast horizon (required when \code{do_forecast = TRUE}).
#' @param X_all    nxpx(T+tnew) covariate array for forecasting (required when
#'   \code{do_forecast = TRUE}).
#' @param spatial  List for spatial interpolation:
#'   \code{list(crd, crdtilde, Xtilde, t)}.
#'   \code{Xtilde} should be a uxp matrix (covariates at prediction sites).
#' @param n_threads  Number of OpenMP threads. Default 1.
#' @param verbose    Logical; print progress? Default TRUE.
#'
#' @return A list with (depending on flags):
#'   \describe{
#'     \item{\code{FF}}{Forward filter output from \code{spFF3_v2}.}
#'     \item{\code{Wi}}{nxJ dynamic weight matrix.}
#'     \item{\code{Wglobal}}{Jx1 global (averaged) weight vector.}
#'     \item{\code{BS}}{(if \code{do_BS}) T-list of nxqxL sample cubes.}
#'     \item{\code{forecast}}{(if \code{do_forecast}) List with \code{Y_pred}
#'       array of dimensions nxqxLxT.}
#'     \item{\code{spatial}}{(if \code{do_spatial}) Named list of L-lists of
#'       prediction matrices, one per requested time.}
#'   }
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 50; T <- 10; p <- 2; q <- 2
#'
#' coords <- matrix(runif(n * 2), ncol = 2)
#' D      <- as.matrix(dist(coords))
#' Y      <- array(rnorm(n * q * T), dim = c(n, q, T))
#' X      <- array(rnorm(n * p * T), dim = c(n, p, T))
#'
#' prior  <- make_prior(
#'   m0   = matrix(0, n + p, q),
#'   C0   = diag(n + p),
#'   nu0  = 3,
#'   Psi0 = diag(q),
#'   p    = p
#' )
#'
#' res <- spFFBS(
#'   Y = Y, X = X, D = D,
#'   grid  = list(tau = c(0.8, 0.9), phi = c(1, 2)),
#'   prior = prior,
#'   L = 50
#' )
#' }
#'
#' @export
spFFBS <- function(
    Y,
    X,
    D,
    grid  = list(tau = NULL, phi = NULL),
    prior,
    G_beta    = NULL,
    rho       = 1.0,
    do_BS     = FALSE,
    do_forecast = FALSE,
    do_spatial  = FALSE,
    L         = 200L,
    tnew      = NULL,
    X_all     = NULL,
    spatial   = NULL,
    n_threads = 1L,
    verbose   = TRUE
) {

  # -- Banner -----------------------------------------------------------------
  if (verbose) {
    cat("\n====================================================\n")
    cat("      Welcome to spFFBS Bayesian Engine\n")
    cat("====================================================\n\n")
  }

  # -- Validate grid ----------------------------------------------------------
  if (is.null(grid$tau) || is.null(grid$phi))
    stop("'grid' must be a list(tau = ..., phi = ...)")

  tau_v <- as.numeric(grid$tau)
  phi_v <- as.numeric(grid$phi)

  if (any(is.na(tau_v)) || any(tau_v <= 0) || any(tau_v >= 1))
    stop(sprintf(
      paste0("All grid$tau values must be in the open interval (0, 1).\n",
             "  Received: %s\n",
             "  Tip: check for typos (e.g. 9975 instead of 0.9975)."),
      paste(tau_v, collapse = ", ")))

  if (any(is.na(phi_v)) || any(phi_v <= 0))
    stop(sprintf(
      "All grid$phi values must be strictly positive. Received: %s",
      paste(phi_v, collapse = ", ")))

  par_grid <- spBPS:::expand_grid_cpp(rev(tau_v), rev(phi_v))
  J    <- nrow(par_grid)
  tmax <- dim(Y)[3L]
  n    <- dim(Y)[1L]
  p_dim <- dim(X)[2L]

  if (verbose)
    cat(sprintf("Parameter grid: %d models | n=%d | T=%d | p=%d\n\n",
                J, n, tmax, p_dim))

  # -- Resolve G_beta ---------------------------------------------------------
  if (is.null(G_beta)) {
    G_beta <- diag(p_dim)   # Identity - most common case
  } else {
    if (!is.matrix(G_beta) || nrow(G_beta) != p_dim || ncol(G_beta) != p_dim)
      stop(sprintf("G_beta must be a %d x %d matrix.", p_dim, p_dim))
  }

  # -- Resolve prior (auto-convert full C0 = blocks) -------------------------
  if (!is.null(prior$C) && is.null(prior$C_bb)) {
    if (verbose) cat("  Auto-converting full prior C to block form...\n")
    prior <- make_prior(prior$m, prior$C, prior$nu, prior$Psi, p_dim)
  }
  required_fields <- c("m", "C_bb", "C_bw", "C_ww", "nu", "Psi")
  missing_f <- setdiff(required_fields, names(prior))
  if (length(missing_f))
    stop("prior is missing fields: ", paste(missing_f, collapse = ", "),
         ". Use make_prior() to build it.")

  # -- Propagate R seed into C++ RNG pool ------------------------------------
  .seed_cpp_rng()

  # -- Forward Filtering ------------------------------------------------------
  if (verbose) cat("[ 1/4 ] Forward Filtering ...\n")
  t0 <- proc.time()["elapsed"]

  out_FF <- spFF3_v2(
    Y        = Y,
    X        = X,
    D        = D,
    par_grid = par_grid,
    prior    = prior,
    G_beta   = G_beta,
    rho      = rho,
    num_threads = as.integer(n_threads)
  )

  dt_FF <- proc.time()["elapsed"] - t0
  if (verbose) cat(sprintf("       Done (%.2f s)\n\n", dt_FF))

  # -- Stacking Weights -------------------------------------------------------
  if (verbose) cat("[ 2/4 ] Computing stacking weights ...\n")
  t0 <- proc.time()["elapsed"]

  Wi      <- compute_Wt(out_FF, tau = grid$tau, phi = grid$phi,
                           n_threads = n_threads, verbose = FALSE)
  Wglobal <- matrix(colMeans(Wi))

  dt_Wi <- proc.time()["elapsed"] - t0
  if (verbose) cat(sprintf("       Done (%.2f s)\n\n", dt_Wi))

  # -- Assemble base result ---------------------------------------------------
  result <- list(FF = out_FF, Wi = Wi, Wglobal = Wglobal)

  # ===========================================================================
  # OPTIONAL: Backward Sampling
  # ===========================================================================
  if (do_BS) {
    if (verbose) cat("[ 3/4 ] Backward Sampling ...\n")
    t0 <- proc.time()["elapsed"]

    out_BS <- weighted_backward_sample_T_v2(
      G_beta    = G_beta,
      rho       = rho,
      D         = D,
      ForwFilt  = out_FF,
      L         = as.integer(L),
      par_grid  = par_grid,
      weights   = as.numeric(Wglobal),
      num_threads = as.integer(n_threads)
    )

    dt_BS <- proc.time()["elapsed"] - t0
    if (verbose) cat(sprintf("       Done (%.2f s)\n\n", dt_BS))
    result$BS <- out_BS
  } else {
    if (verbose) cat("[ 3/4 ] Backward Sampling skipped.\n\n")
  }

  # ===========================================================================
  # OPTIONAL: Temporal Forecast
  # ===========================================================================
  if (do_forecast) {
    if (is.null(tnew) || tnew <= 0L)
      stop("do_forecast = TRUE requires tnew > 0.")
    if (is.null(X_all))
      stop("do_forecast = TRUE requires X_all: the nxpx(T+tnew) covariate array.")
    if (dim(X_all)[3L] < tmax + tnew)
      stop(sprintf("X_all must have at least %d time slices (T + tnew = %d + %d).",
                   tmax + tnew, tmax, tnew))

    if (verbose) cat("[ 4/4 ] Temporal Forecast ...\n")
    t0 <- proc.time()["elapsed"]

    # Dynamic weight matrix: J x (tmax-1), each col = Wglobal
    Wglobal_dyn <- matrix(rep(as.numeric(Wglobal), tmax - 1L), nrow = J)

    out_TF <- temporal_forecast_v2(
      G_beta    = G_beta,
      rho       = rho,
      D         = D,
      par_grid  = par_grid,
      ForwFilt  = out_FF,
      X_all     = X_all,
      weights   = Wglobal_dyn,
      horiz     = as.integer(tnew),
      L         = as.integer(L),
      num_threads = as.integer(n_threads)
    )

    Y_pred  <- abind::abind(out_TF, along = 4L)
    dt_TF   <- proc.time()["elapsed"] - t0
    if (verbose) cat(sprintf("       Done (%.2f s)\n\n", dt_TF))
    result$forecast <- list(out_TF = out_TF, Y_pred = Y_pred)
  } else {
    if (verbose) cat("[ 4/4 ] Temporal Forecast skipped.\n\n")
  }

  # ===========================================================================
  # OPTIONAL: Spatial Interpolation
  # ===========================================================================
  if (do_spatial) {
    if (is.null(spatial))
      stop("do_spatial = TRUE requires spatial = list(crd, crdtilde, Xtilde, t).")

    required_sp <- c("crd", "crdtilde", "Xtilde", "t")
    missing_sp  <- setdiff(required_sp, names(spatial))
    if (length(missing_sp))
      stop("spatial is missing: ", paste(missing_sp, collapse = ", "))

    crd      <- spatial$crd
    crdtilde <- spatial$crdtilde
    Xtilde   <- spatial$Xtilde   # u x p matrix OR u x p x length(t) array
    t_vec    <- as.integer(spatial$t)

    u <- nrow(crdtilde)

    # -- Validate and normalise Xtilde -------------------------------------
    # Accept either:
    #   (a) u x p matrix  - same covariates used for every t in t_vec
    #   (b) u x p x K array - one covariate slice per element of t_vec
    if (is.array(Xtilde) && length(dim(Xtilde)) == 3L) {
      if (dim(Xtilde)[1L] != u || dim(Xtilde)[2L] != p_dim)
        stop(sprintf(
          "spatial$Xtilde array must have dimensions u=%d x p=%d x K. Got %s.",
          u, p_dim, paste(dim(Xtilde), collapse = " x ")))
      K_sp <- dim(Xtilde)[3L]
      if (K_sp != length(t_vec))
        stop(sprintf(
          paste0("spatial$Xtilde has %d slices but spatial$t has %d elements.\n",
                 "  They must match: one covariate slice per time point."),
          K_sp, length(t_vec)))
      Xtilde_is_array <- TRUE
    } else {
      if (!is.matrix(Xtilde) || nrow(Xtilde) != u || ncol(Xtilde) != p_dim)
        stop(sprintf(
          paste0("spatial$Xtilde must be a u=%d x p=%d matrix, ",
                 "or a u x p x length(t) array."),
          u, p_dim))
      Xtilde_is_array <- FALSE
    }

    D_us  <- spBPS:::arma_dist(rbind(crdtilde, crd))[seq_len(u), -seq_len(u)]
    D_u   <- spBPS:::arma_dist(crdtilde)

    spatial_out <- list()

    if (verbose) cat("[+]   Spatial Interpolation ...\n")
    t0_sp <- proc.time()["elapsed"]

    for (idx_t in seq_along(t_vec)) {
      t_pick <- t_vec[idx_t]

      if (t_pick < 1L)
        stop(sprintf("spatial$t contains invalid value %d: must be >= 1.", t_pick))

      # Select covariate matrix for this t
      Xu_t <- if (Xtilde_is_array) Xtilde[, , idx_t, drop = FALSE][,,1] else Xtilde

      # Dynamic weight matrix: J x min(t_pick, tmax) columns.
      t_w <- min(t_pick, tmax)
      W_dyn <- cbind(
        matrix(rep(1.0 / J, J), ncol = 1L),
        matrix(rep(as.numeric(Wglobal), t_w - 1L), nrow = J)
      )

      out_SI <- spatial_interpolation_v2(
        G_beta    = G_beta,
        rho       = rho,
        Xu        = Xu_t,
        D_s       = D,
        D_u       = D_u,
        D_us      = D_us,
        par_grid  = par_grid,
        ForwFilt  = out_FF,
        weights   = W_dyn,
        t         = as.integer(t_pick),
        L         = as.integer(L),
        num_threads = as.integer(n_threads)
      )

      # Convert list of L matrices (each (2u)xq) to two uxqxL arrays.
      # vapply pre-allocates the result buffer - no repeated reallocation.
      q_dim <- ncol(out_SI[[1L]])
      sp_arr <- vapply(out_SI, as.numeric, numeric(2L * u * q_dim))
      # sp_arr is (2u*q) x L column-major; reshape to (2u, q, L)
      sp_arr <- array(sp_arr, dim = c(2L * u, q_dim, L))
      spatial_out[[paste0("t", t_pick)]] <- list(
        Y     = sp_arr[seq_len(u),     , , drop = FALSE],
        Omega = sp_arr[u + seq_len(u), , , drop = FALSE]
      )
    }

    dt_sp <- proc.time()["elapsed"] - t0_sp
    if (verbose) cat(sprintf("       Done (%.2f s)\n\n", dt_sp))
    result$spatial <- spatial_out
  }

  # -- Summary ----------------------------------------------------------------
  if (verbose) {
    cat("====================================================\n")
    cat("    spFFBS pipeline completed successfully!\n")
    cat("====================================================\n\n")
  }

  result
}
