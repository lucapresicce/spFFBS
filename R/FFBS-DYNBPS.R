

#' compute_Wt: Dynamic Bayesian Predictive Stacking Weights
#'
#' Computes predictive stacking weights from the forward-filter objects.
#' Provides a simplified and user-friendly interface with optional parallelization.
#'
#' @param out_FF Output of spFF3 (list of filtering results).
#' @param tau Vector of tau grid values (only used for column naming).
#' @param phi Vector of phi grid values (only used for column naming).
#' @param parallel Logical; use parallel backend (foreach + doParallel)? Default FALSE.
#'
#' @return Matrix of weights of size n x J
#' @keywords internal
compute_Wt <- function(out_FF, tau, phi, parallel = FALSE) {

  cat("\n====================================================\n")
  cat("    ⚖️  Computing Dynamic Bayesian Predictive\n")
  cat("             Stacking Weights (Wi)\n")
  cat("====================================================\n\n")

  tictoc::tic("Total time")

  # ---------------------------------------------------------------------------
  # Extract basic filter info
  # ---------------------------------------------------------------------------
  tmax <- length(out_FF)
  if (tmax < 2) stop("Need at least two time slices for predictive stacking.")

  density0 <- out_FF[[2]]$density_evaluations
  if (is.null(density0))
    stop("Forward filter object does not contain $density_evaluations.")

  n  <- nrow(density0)
  J  <- ncol(density0)

  cat("📦 Found:\n")
  cat("   • Time points:  ", tmax, "\n")
  cat("   • Locations:    ", n, "\n")
  cat("   • Models (J):   ", J, "\n\n")

  # ---------------------------------------------------------------------------
  # Choose backend
  # ---------------------------------------------------------------------------
  if (parallel) {
    cat("🔧 Using PARALLEL backend (foreach + registered cluster)\n")
    `%run%` <- foreach::`%dopar%`
  } else {
    cat("🔧 Using SERIAL backend (foreach sequential)\n")
    `%run%` <- foreach::`%do%`
  }
  cat("\n")

  # ---------------------------------------------------------------------------
  # MAIN LOOP over locations
  # ---------------------------------------------------------------------------
  cat("⚙️  Evaluating model scores and optimizing weights...\n")
  tictoc::tic("Weight computation")

  Wi <- foreach::foreach(i = 1:n, .combine = "rbind", .packages = "spBPS") %run% {

    # Collect log predictive densities over times for location i
    epd_i <- sapply(2:tmax, function(a) out_FF[[a]]$density_evaluations[i, ])

    # Handle degenerate cases
    if (is.null(dim(epd_i))) {
      # No meaningful variation → uniform weights
      w <- rep(1 / J, J)
      return(matrix(w, nrow = 1))
    }

    # Scores for stacking (exponentiate log predictive densities)
    scores <- exp(t(epd_i))   # J x (tmax-1)

    # Robust fallback for zeros or degenerate scores
    if (all(scores == 0) || anyNA(scores)) {
      w <- rep(1 / J, J)
      return(matrix(w, nrow = 1))
    }

    # Optimization (your internal function)
    Wi_i <- try(optimize_weights_proj(scores = scores), silent = TRUE)

    if (inherits(Wi_i, "try-error")) {
      # Safe fallback
      w <- rep(1 / J, J)
      return(matrix(w, nrow = 1))
    }

    # Return row vector
    matrix(as.numeric(Wi_i), nrow = 1)
  }

  tictoc::toc()

  # ---------------------------------------------------------------------------
  # Name the models
  # ---------------------------------------------------------------------------
  colnames(Wi) <- paste0(
    "Model ", 1:J,
    "  [ tau=", rev(tau), ", phi=", rev(phi), " ]"
  )

  cat("\n✔️  Weight matrix computed successfully.\n")
  cat("   Dimensions: ", nrow(Wi), "x", ncol(Wi), "\n\n")
  cat("====================================================\n\n")

  tictoc::toc()

  return(Wi)
}



#' spFFBS: Spatiotemporal Bayesian Pipeline (friendly interface)
#'
#' A user-friendly, modular wrapper for running a Bayesian spatiotemporal
#' filtering + weighting pipeline, with optional backward sampling,
#' forecasting, and spatial interpolation.
#'
#' @param Y Response data (3D array or cube).
#' @param G System matrix (cube).
#' @param P Observation matrix (cube).
#' @param D Spatial distance matrix.
#' @param grid List with elements:
#'   - tau: numeric vector
#'   - phi: numeric vector
#'
#' @param prior Prior list for forward filter (m, C, nu, Psi)
#'
#' @param do_BS Logical: run backward sampling? (default: FALSE)
#' @param do_forecast Logical: run temporal forecasts? (default: FALSE)
#' @param do_spatial Logical: run spatial interpolation? (default: FALSE)
#'
#' @param L Number of posterior samples (default 200)
#' @param tnew Forecast horizon (default 5)
#'
#' @param spatial Optional list for spatial:
#'   list(crd = , crdtilde = , Xtilde = )
#'
#' @return A list with the components executed according to the flags.
#' @export
spFFBS <- function(
    Y, G, P, D,
    grid = list(tau = NULL, phi = NULL),
    prior,
    do_BS = FALSE,
    do_forecast = FALSE,
    do_spatial = FALSE,
    L = 200,
    tnew = NULL,
    spatial = NULL
) {

  cat("\n====================================================\n")
  cat("        🚀 Welcome to spFFBS Bayesian Engine\n")
  cat("====================================================\n\n")

  # --- Parameter grid --------------------------------------------------------
  if (is.null(grid$tau) || is.null(grid$phi)) {
    stop("'grid' must be a list(tau = ..., phi = ...)")
  }

  cat("🔧 Building parameter grid ... ")
  par_grid <- spBPS::expand_grid_cpp(rev(grid$tau), rev(grid$phi))
  J <- nrow(par_grid)
  cat("OK (", J, "models )\n\n")

  tmax <- dim(Y)[3]
  n <- dim(Y)[1]

  # --- Forward Filtering -----------------------------------------------------
  cat("📡 Running Forward Filtering (FF)...\n")
  tictoc::tic()
  out_FF <- spFF3(
    Y = Y, G = G, P = P, D = D,
    par_grid = par_grid, prior = prior, num_threads = 1
  )
  tictoc::toc()
  cat("✔️ FF completed.\n\n")

  # --- Weight Computation ----------------------------------------------------
  cat("⚖️ Computing stacking weights ...\n")
  tictoc::tic()
  Wi <- compute_Wt(out_FF, tau = grid$tau, phi = grid$phi)
  tictoc::toc()
  Wglobal <- matrix(colMeans(Wi))
  cat("✔️ Global weights computed.\n\n")

  # Initialize result list
  result <- list(
    FF = out_FF,
    Wi = Wi,
    Wglobal = Wglobal
  )

  # ---------------------------------------------------------------------------
  # OPTIONAL PARTS
  # ---------------------------------------------------------------------------

  # === Backward Sampling =====================================================
  if (do_BS) {
    cat("🎯 Running Backward Sampling (BS)...\n")
    tictoc::tic()

    repeat {
      res_try <- try({
        out_BS <- weighted_backward_sample_T3(
          G = G, D = D, ForwFilt = out_FF,
          L = L, par_grid = par_grid, weights = Wglobal
        )
      }, silent = TRUE)

      if (!inherits(res_try, "try-error")) break
      message("   ⚠️  Retry due to numerical error...")
    }

    tictoc::toc()
    cat("✔️ Backward sampling completed.\n\n")
    result$BS <- out_BS
  }

  # === Forecasting ===========================================================
  if (do_forecast) {
    if (is.null(tnew)) stop("If (do_forecast = TRUE), you must supply tnew > 0.")
    if (dim(G)[3] != tmax + tnew || dim(P)[3] != tmax + tnew) {
      stop("For forecasting (do_forecast = TRUE),
           G and P must contain (tmax + tnew) slices covering the out-of-sample period.")
    }

    cat("🔮 Running Temporal Forecast (TF)...\n")
    Wglobal_dyn <- matrix(rep(as.numeric(Wglobal), tmax - 1), nrow = J)

    tictoc::tic()
    out_TF <- temporal_forecast3(
      G = G, P = P, D = D,
      par_grid = par_grid,
      ForwFilt = out_FF,
      weights = Wglobal_dyn,
      horiz = tnew, L = L
    )
    tictoc::toc()
    Y_pred <- abind::abind(out_TF, along = 4)
    cat("✔️ Forecasting completed.\n\n")

    result$forecast <- list(out_TF = out_TF, Y_pred = Y_pred)
  }

  # === Spatial Interpolation =================================================
  if (do_spatial) {
    if (is.null(spatial)) {
      stop("If (do_spatial = TRUE), provide spatial=list(crd=, crdtilde=, Xtilde=, t=).")
    }

    # Validate time index
    if (is.null(spatial$t)) {
      stop("When (do_spatial = TRUE), you must specify spatial$t (time index).")
    }

    t_vec <- spatial$t
    if (!is.numeric(t_vec)) stop("spatial$t must be numeric.")
    if (any(t_vec < 1) || any(t_vec > dim(G)[3])) {
      stop("spatial$t contains times outside available G slices.")
    }

    crd      <- spatial$crd
    crdtilde <- spatial$crdtilde
    Xtilde   <- spatial$Xtilde
    u <- nrow(crdtilde)

    # Distances
    D_us <- spBPS::arma_dist(rbind(crdtilde, crd))[1:u, -(1:u)]
    D_tilde <- spBPS::arma_dist(crdtilde)

    # Prepare outputs
    spatial_out <- list()

    cat("🗺️ Running Spatial Interpolation ...\n")

    tictoc::tic()
    for (t_pick in t_vec) {

      # dynamic weight matrix
      W_dyn <- cbind(
        matrix(rep(1 / J, J), ncol = 1),
        matrix(rep(as.numeric(Wglobal), t_pick - 1), nrow = J)
      )

      # replicate last FF if needed
      FF_t <- out_FF
      if (t_pick > length(out_FF)) {
        FF_t[[t_pick]] <- out_FF[[tmax]]
      }

      out_SI <- spatial_interpolation(
        G = G[,,t_pick], P = P[,,t_pick], Xu = Xtilde[,,t_pick],
        D_s = D, D_u = D_tilde, D_us = D_us,
        par_grid = par_grid, ForwFilt = FF_t,
        weights = W_dyn, t = t_pick, L = L
      )

      spatial_out[[ paste0("t", t_pick) ]] <- out_SI
    }
    tictoc::toc()

    result$spatial <- spatial_out
    cat("✔️ Spatial interpolation completed.\n\n")
  }

  cat("====================================================\n")
  cat("     🎉 spFFBS pipeline completed successfully!\n")
  cat("====================================================\n\n")

  return(result)
}
