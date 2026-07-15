// ============================================================
//  FFBS-DYNBPS-struct-v2.cpp  (revision 2)
//  Optimised C++ backend for the spFFBS package.
// ============================================================

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <atomic>
#include <random>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;


// ============================================================
// SECTION 1 - Thread-local RNG pool
// ============================================================
namespace {

struct EnginePool_v2 {
  std::atomic<uint64_t> base_seed{5489u};
  std::atomic<uint64_t> counter{0u};
  std::atomic<uint64_t> epoch{0u};

  static EnginePool_v2& instance() {
    static EnginePool_v2 pool;
    return pool;
  }
  void reseed(uint64_t seed) {
    if (seed == 0u) seed = 1u;
    base_seed.store(seed, std::memory_order_relaxed);
    counter.store(0u, std::memory_order_relaxed);
    epoch.fetch_add(1u, std::memory_order_acq_rel);
  }
  uint64_t next_seed() {
    uint64_t id = counter.fetch_add(1u, std::memory_order_relaxed) + 1u;
    constexpr uint64_t phi = 0x9E3779B97F4A7C15ULL;
    return base_seed.load(std::memory_order_relaxed) + phi * id;
  }
  uint64_t current_epoch() const {
    return epoch.load(std::memory_order_acquire);
  }
};

inline std::mt19937_64& thread_engine() {
  // Seed deterministically by thread ID, not by atomic counter order.
  // This guarantees that with the same base seed and the same n_threads,
  // results are identical across runs - even if threads initialize in
  // different orders.  Different n_threads - different per-thread seeds
  // - different (but still valid) Monte Carlo draws: document in Rd.
  thread_local std::mt19937_64 eng(0u);
  thread_local uint64_t eng_epoch = std::numeric_limits<uint64_t>::max();
  uint64_t pool_epoch = EnginePool_v2::instance().current_epoch();
  if (eng_epoch != pool_epoch) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    // Each thread gets a unique, deterministic seed:
    //   base_seed  XOR  phi_hash(tid+1)
    // phi = 2^64 / golden_ratio - maximally separated seeds in [0, 2^64)
    constexpr uint64_t phi = 0x9E3779B97F4A7C15ULL;
    uint64_t ts = EnginePool_v2::instance().base_seed.load(
      std::memory_order_relaxed)
      ^ (phi * static_cast<uint64_t>(tid + 1));
    eng.seed(ts);
    eng_epoch = pool_epoch;
  }
  return eng;
}

inline std::normal_distribution<double>& normal_dist() {
  thread_local std::normal_distribution<double> d(0.0, 1.0);
  return d;
}

inline arma::mat std_normal_mat(arma::uword nr, arma::uword nc) {
  arma::mat Z(nr, nc);
  auto& eng = thread_engine();
  auto& d   = normal_dist();
  for (arma::uword j = 0; j < nc; ++j)
    for (arma::uword i = 0; i < nr; ++i)
      Z(i, j) = d(eng);
  return Z;
}

inline double std_normal_scalar() { return normal_dist()(thread_engine()); }

inline double chi2_scalar(double df) {
  std::gamma_distribution<double> gd(0.5 * df, 2.0);
  return gd(thread_engine());
}

} // anonymous namespace


// [[Rcpp::export]]
void seed_rng_v2(double u1, double u2) {
  uint64_t s1 = static_cast<uint64_t>(u1 * 4294967295.0);
  uint64_t s2 = static_cast<uint64_t>(u2 * 4294967295.0);
  uint64_t seed = (s1 << 32) ^ s2 ^ 0x9E3779B97F4A7C15ULL;
  EnginePool_v2::instance().reseed(seed);
}


// ============================================================
// SECTION 2 - Utilities
// ============================================================

template<class Fn>
static void spFFBS_par_for(int K, int n_cores, Fn fn) {
  if (n_cores <= 1) {
    for (int k = 0; k < K; ++k) fn(k);
    return;
  }
#ifdef _OPENMP
  omp_set_max_active_levels(1);
#pragma omp parallel for schedule(static) num_threads(n_cores)
  for (int k = 0; k < K; ++k) fn(k);
#else
  for (int k = 0; k < K; ++k) fn(k);
#endif
}

// safe_chol: jitter fallback - kept only for posterior covariances.
// HARDENED: if all Cholesky attempts fail (e.g. due to invalid tau > 1
// causing V < 0), return a diagonal fallback rather than an indeterminate
// matrix.  All callers check Q_diag <= 0 and set density = -1e15, so a
// fallback L_Q = sqrt(diag(|S|)) still produces valid (if wide) samples
// and never causes a dimension-mismatch abort.
static arma::mat safe_chol(const arma::mat& M) {
  int n = static_cast<int>(M.n_rows);
  arma::mat S = 0.5 * (M + M.t());
  arma::mat L;
  if (arma::chol(L, S, "lower")) return L;
  double jit = 1e-8 * arma::trace(S) / n;
  if (jit <= 0.0) jit = 1e-8;
  S.diag() += jit;
  if (arma::chol(L, S, "lower")) return L;
  S.diag() += 1e-4 * arma::norm(M, "fro") / n;
  if (arma::chol(L, S, "lower")) return L;
  // Last-resort fallback: diagonal matrix from |diag(M)|.
  // This keeps L non-empty (correct size) so no subsequent
  // matrix-product dimension mismatch can cause an abort.
  L.zeros(n, n);
  for (int i = 0; i < n; ++i)
    L(i, i) = std::sqrt(std::max(std::abs(M(i, i)), 1e-10));
  return L;
}

// -- Pure-C++ Cholesky and SPD inverse ------------------------
// These functions call NO BLAS/LAPACK internally, so they are safe to
// call from concurrent OMP threads regardless of which BLAS vendor
// (OpenBLAS, MKL, BLIS) is installed.  OpenBLAS switches to its internal
// thread pool for matrices >= ~128x128; calling it from 5 OMP threads
// causes over-subscription and deterministic abort().  These routines
// replace arma::inv_sympd / safe_chol for the n x n and np x np
// operations inside the parallel cache-building loop.
// p x p operations (p <= 8) keep arma::inv_sympd - BLAS never spawns
// threads for 8x8 matrices.

// Lower Cholesky: fills out with L such that L * L^T = A.
// Upper triangle of out is zeroed.  Returns false if A is not SPD.
static bool chol_lower_noblas(arma::mat& out, const arma::mat& A) {
  const int n = static_cast<int>(A.n_rows);
  out.set_size(n, n);         // no-op if same size (pre-allocated)
  for (int j = 0; j < n; ++j) {
    // Zero upper triangle of column j
    for (int i = 0; i < j; ++i) out(i, j) = 0.0;
    // Diagonal
    double s = A(j, j);
    for (int k = 0; k < j; ++k) s -= out(j, k) * out(j, k);
    if (s <= 0.0) return false;
    out(j, j) = std::sqrt(s);
    const double inv_d = 1.0 / out(j, j);
    // Sub-diagonal entries of column j
    for (int i = j + 1; i < n; ++i) {
      double t = A(i, j);
      for (int k = 0; k < j; ++k) t -= out(i, k) * out(j, k);
      out(i, j) = t * inv_d;
    }
  }
  return true;
}

// Invert SPD matrix A.  Algorithm:
//   1. L = chol(A)             [lower triangular, pure C++]
//   2. Li = inv(L)             [forward substitution, O(n^3/3)]
//   3. A^{-1} = Li^T * Li     [symmetric rank-n, O(n^3/2)]
// Total: O(n^3) with constant ~1 - comparable to LAPACK potrf+potri.
// Allocates two local n x n matrices (Li and workspace in L); these are
// per-thread stack objects and do not cause malloc contention.
[[maybe_unused]] static bool inv_sympd_noblas(arma::mat& out, const arma::mat& A) {
  const int n = static_cast<int>(A.n_rows);
  // Step 1: Cholesky into local L
  arma::mat L;
  if (!chol_lower_noblas(L, A)) return false;
  // Step 2: Li = L^{-1} via forward substitution (Li is lower triangular)
  //   Li(j,j) = 1/L(j,j)
  //   Li(i,j) = -[ sum_{k=j}^{i-1} L(i,k)*Li(k,j) ] / L(i,i)   i > j
  arma::mat Li(n, n, arma::fill::zeros);
  for (int j = 0; j < n; ++j) {
    Li(j, j) = 1.0 / L(j, j);
    for (int i = j + 1; i < n; ++i) {
      double s = 0.0;
      for (int k = j; k < i; ++k) s += L(i, k) * Li(k, j);
      Li(i, j) = -s / L(i, i);
    }
  }
  // Step 3: A^{-1}(i,j) = sum_{k=max(i,j)}^{n-1} Li(k,i)*Li(k,j)
  // Li is lower triangular, so Li(k,i) = 0 for k < i.
  out.set_size(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double s = 0.0;
      for (int k = i; k < n; ++k) s += Li(k, i) * Li(k, j);
      out(i, j) = out(j, i) = s;
    }
  }
  return true;
}
// -- End pure-C++ helpers -----------------------------------------------------

static double log_mvgamma(double alpha, int q) {
  double v = 0.25 * static_cast<double>(q * (q - 1)) * std::log(M_PI);
  for (int i = 0; i < q; ++i)
    v += std::lgamma(alpha - 0.5 * static_cast<double>(i));
  return v;
}


// ============================================================
// SECTION 3 - C++ samplers
// ============================================================

static arma::mat rIW(double nu, const arma::mat& Psi) {
  arma::uword q   = Psi.n_rows;
  arma::mat L_psi = safe_chol(Psi);
  arma::mat L_inv = arma::solve(arma::trimatl(L_psi),
                                arma::eye<arma::mat>(q, q));
  arma::mat A(q, q, arma::fill::zeros);
  for (arma::uword i = 0; i < q; ++i) {
    A(i, i) = std::sqrt(chi2_scalar(nu - static_cast<double>(i)));
    for (arma::uword j = i + 1; j < q; ++j)
      A(j, i) = std_normal_scalar();
  }
  arma::mat L = L_inv * A;
  arma::mat W = 0.5 * (L * L.t() + (L * L.t()).t());
  arma::mat X;
  if (!arma::inv_sympd(X, W)) {
    W.diag() += 1e-10 * arma::trace(W) / q;
    arma::inv_sympd(X, W);
  }
  return X;
}

// L_C (chol of C) passed in - no re-factorisation per draw
static void rMNIW_batch_Lc(
    arma::uword n_draws,
    const arma::mat& m,
    const arma::mat& L_C,
    double nu, const arma::mat& Psi,
    arma::cube& Theta_out,
    arma::cube& Sigma_out,
    int n_threads = 1) {
  arma::uword nr = m.n_rows, nc = m.n_cols;
  Theta_out.set_size(nr, nc, n_draws);
  Sigma_out.set_size(nc, nc, n_draws);
  spFFBS_par_for(static_cast<int>(n_draws), n_threads,
                 [&](int l) {
                   arma::mat Sig   = rIW(nu, Psi);
                   arma::mat L_Sig = safe_chol(Sig);
                   arma::mat Z     = std_normal_mat(nr, nc);
                   Theta_out.slice(l) = m + L_C * Z * L_Sig.t();
                   Sigma_out.slice(l) = Sig;
                 });
}


// ============================================================
// SECTION 4 - Matrix-T log-density
// ============================================================

// [[Rcpp::export]]
double log_matrix_t_density_v2(
    const arma::mat& X,
    const arma::mat& M,
    double U_scalar,
    const arma::mat& L_V,
    double log_det_V,
    double lgamma_nupq,
    double lgamma_nuq,
    double nu) {
  int q = X.n_cols;
  double log_det_U = std::log(U_scalar);
  arma::mat diff   = X - M;
  arma::mat A      = arma::solve(arma::trimatl(L_V), diff.t());
  double B_scalar  = U_scalar + arma::as_scalar(A.t() * A);
  double ldet_ratio = 0.5 * (std::log(B_scalar) - log_det_U);
  return lgamma_nupq - lgamma_nuq
  - (nu + static_cast<double>(q)) * ldet_ratio
  - static_cast<double>(q) * 0.5 * log_det_U
  - 0.5 * log_det_V
  - static_cast<double>(q) * 0.5 * std::log(M_PI);
}


// ============================================================
// SECTION 5 - FFResult_v2
// ============================================================

struct FFResult_v2 {
  // Backward sampling
  arma::mat m;        // (p+n) x q  filtered mean
  arma::mat C_bb;     // p x p
  arma::mat C_bw;     // p x n
  arma::mat C_ww;     // n x n  (dominant)
  arma::mat Psi;      // q x q
  double    nu;
  // Density collapse
  arma::vec Q_diag;   // n x 1  exact diagonal of Q
  arma::mat f;        // n x q  predictive mean
  arma::mat L_Psi;    // q x q  lower Chol of Psi (precomputed)
  double log_det_Psi;
  double lgmm_nupq;
  double lgmm_nuq;
};


// ============================================================
// SECTION 6 - Block-structured forward filter
// ============================================================

static FFResult_v2 forward_filter_v2(
    const arma::mat& Y,
    const arma::mat& X,
    const arma::mat& G_beta,
    double rho,
    bool   G_is_identity,
    const arma::mat& K_j,    // precomputed
    double c_j,
    const arma::mat& prev_m,
    const arma::mat& prev_C_bb,
    const arma::mat& prev_C_bw,
    const arma::mat& prev_C_ww,
    double nu,
    const arma::mat& Psi) {

  int n = Y.n_rows, q = Y.n_cols, p = G_beta.n_rows;

  arma::mat m_b = prev_m.rows(0, p - 1);
  arma::mat m_w = prev_m.rows(p, p + n - 1);

  // -- Prediction step ---------------------------------------------
  arma::mat a_b, a_w, R_bb, R_bw, R_ww;
  if (G_is_identity) {
    a_b = m_b;  a_w = m_w;
    R_bb = prev_C_bb;  R_bb.diag() += 1.0;
    R_bw = prev_C_bw;
    R_ww = prev_C_ww + K_j;
  } else {
    a_b  = G_beta * m_b;
    a_w  = rho * m_w;
    R_bb = G_beta * prev_C_bb * G_beta.t();  R_bb.diag() += 1.0;
    R_bw = rho * (G_beta * prev_C_bw);
    R_ww = (rho * rho) * prev_C_ww + K_j;
  }

  // -- Q = A_bot + X*A_top + c*I  ----------
  arma::mat A_top = R_bb * X.t() + R_bw;       // p x n
  arma::mat A_bot = R_bw.t() * X.t() + R_ww;   // n x n

  arma::mat Q = A_bot + X * A_top;
  Q.diag() += c_j;
  Q = 0.5 * (Q + Q.t());

  // Cholesky - Q is strictly SPD by construction
  arma::mat L_Q;
  if (!arma::chol(L_Q, Q, "lower")) L_Q = safe_chol(Q);

  // -- Innovation and solves ----------------------------------------
  arma::mat f_pred = X * a_b + a_w;
  arma::mat e      = Y - f_pred;

  arma::mat W_e_fwd   = arma::solve(arma::trimatl(L_Q), e);
  arma::mat W_e       = arma::solve(arma::trimatu(L_Q.t()), W_e_fwd);
  arma::mat W_top_fwd = arma::solve(arma::trimatl(L_Q), A_top.t());
  arma::mat W_top     = arma::solve(arma::trimatu(L_Q.t()), W_top_fwd);
  arma::mat W_bot_fwd = arma::solve(arma::trimatl(L_Q), A_bot.t());

  // -- Covariance update ----------------------------------
  FFResult_v2 res;
  arma::mat C_bb_new = R_bb - A_top * W_top;
  res.C_bb = 0.5 * (C_bb_new + C_bb_new.t());
  res.C_bw = R_bw - W_top.t() * A_bot.t();
  arma::mat C_ww_new = R_ww - W_bot_fwd.t() * W_bot_fwd;
  res.C_ww = 0.5 * (C_ww_new + C_ww_new.t());

  // -- Mean and sufficient statistics ------------------------------
  res.m.set_size(p + n, q);
  res.m.rows(0, p - 1)     = a_b + W_top.t() * e;
  res.m.rows(p, p + n - 1) = a_w + A_bot * W_e;

  res.nu  = nu + 0.5 * static_cast<double>(n);
  arma::mat Psi_new = Psi + 0.5 * (e.t() * W_e);
  res.Psi = 0.5 * (Psi_new + Psi_new.t());

  // -- Density fields for collapse(2) ---------------------
  res.Q_diag = Q.diag();
  res.f      = f_pred;
  res.L_Psi  = safe_chol(res.Psi);
  res.log_det_Psi = 2.0 * arma::sum(arma::log(res.L_Psi.diag()));
  double nuq_d   = res.nu + static_cast<double>(q) - 1.0;
  double nupq_d  = nuq_d + 1.0;
  res.lgmm_nupq  = log_mvgamma(0.5 * nupq_d, q);
  res.lgmm_nuq   = log_mvgamma(0.5 * nuq_d,  q);

  return res;
}


// ============================================================
// SECTION 7 - Unified parallel function
// ============================================================

// Internal version used by spFF3_v2 (receives K_vec directly)
static void run_filters_and_density(
    const arma::mat& Y,
    const arma::mat& X,
    const arma::mat& G_beta,
    double rho,
    bool   G_is_identity,
    const std::vector<arma::mat>& K_vec,
    const std::vector<double>&    c_vec,
    const std::vector<arma::mat>& m_vec,
    const std::vector<arma::mat>& C_bb_vec,
    const std::vector<arma::mat>& C_bw_vec,
    const std::vector<arma::mat>& C_ww_vec,
    const std::vector<double>&    nu_vec,
    const std::vector<arma::mat>& Psi_vec,
    int J, int n, int num_threads,
    std::vector<FFResult_v2>& results,
    arma::mat& density_evaluations) {

  results.resize(J);
  density_evaluations.set_size(n, J);
  density_evaluations.zeros();

  // Forward filter - OMP over j
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(num_threads)
#endif
  for (int j = 0; j < J; j++) {
    results[j] = forward_filter_v2(
      Y, X, G_beta, rho, G_is_identity,
      K_vec[j], c_vec[j],
                     m_vec[j], C_bb_vec[j], C_bw_vec[j], C_ww_vec[j],
                                                                 nu_vec[j], Psi_vec[j]);
  }

  // Density - OMP collapse(2) over (j,g)
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) num_threads(num_threads)
#endif
  for (int j = 0; j < J; j++) {
    for (int g = 0; g < n; g++) {
      double Qgg = results[j].Q_diag(g);
      if (Qgg <= 0.0) { density_evaluations(g, j) = -1e15; continue; }
      density_evaluations(g, j) = log_matrix_t_density_v2(
        Y.row(g), results[j].f.row(g), Qgg,
        results[j].L_Psi, results[j].log_det_Psi,
        results[j].lgmm_nupq, results[j].lgmm_nuq,
        results[j].nu);
    }
  }
}


// Exported version (stand-alone R calls) - computes K_vec on-the-fly
// [[Rcpp::export]]
Rcpp::List unified_parallel_function_v2(
    const arma::mat& Y,
    const arma::mat& X,
    const arma::mat& G_beta,
    double rho,
    bool   G_is_identity,
    const arma::mat& D,
    const arma::mat& par_grid,
    Rcpp::List FF_prec,
    int J, int p, int n,
    int num_threads = 1) {

  std::vector<arma::mat> m_vec(J), C_bb_vec(J), C_bw_vec(J), C_ww_vec(J);
  std::vector<arma::mat> Psi_vec(J);
  std::vector<double>    nu_vec(J), c_vec(J);

  for (int j = 0; j < J; j++) {
    Rcpp::List res_j = FF_prec[j];
    m_vec[j]    = Rcpp::as<arma::mat>(res_j["m"]);
    C_bb_vec[j] = Rcpp::as<arma::mat>(res_j["C_bb"]);
    C_bw_vec[j] = Rcpp::as<arma::mat>(res_j["C_bw"]);
    C_ww_vec[j] = Rcpp::as<arma::mat>(res_j["C_ww"]);
    Psi_vec[j]  = Rcpp::as<arma::mat>(res_j["Psi"]);
    nu_vec[j]   = Rcpp::as<double>(res_j["nu"]);
    c_vec[j]    = (1.0 - par_grid(j, 0)) / par_grid(j, 0);
  }

  // Compute K_vec on-the-fly (stand-alone path)
  std::vector<arma::mat> K_vec(J);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(num_threads)
#endif
  for (int j = 0; j < J; j++)
    K_vec[j] = arma::exp(-par_grid(j, 1) * D);

  std::vector<FFResult_v2> results;
  arma::mat density_evaluations;
  run_filters_and_density(Y, X, G_beta, rho, G_is_identity,
                          K_vec, c_vec,
                          m_vec, C_bb_vec, C_bw_vec, C_ww_vec,
                          nu_vec, Psi_vec,
                          J, n, num_threads,
                          results, density_evaluations);

  Rcpp::List filtered_results(J);
  for (int j = 0; j < J; j++) {
    const FFResult_v2& r = results[j];
    filtered_results[j] = Rcpp::List::create(
      Rcpp::Named("m")    = r.m,
      Rcpp::Named("C_bb") = r.C_bb,
      Rcpp::Named("C_bw") = r.C_bw,
      Rcpp::Named("C_ww") = r.C_ww,
      Rcpp::Named("nu")   = r.nu,
      Rcpp::Named("Psi")  = r.Psi);
  }
  return Rcpp::List::create(
    Rcpp::Named("filtered_results")    = filtered_results,
    Rcpp::Named("density_evaluations") = density_evaluations);
}


// ============================================================
// SECTION 8 - Main forward-filtering loop
// ============================================================

// [[Rcpp::export]]
Rcpp::List spFF3_v2(
    const arma::cube& Y,
    const arma::cube& X,
    const arma::mat&  D,
    const arma::mat&  par_grid,
    Rcpp::List const& prior,
    const arma::mat&  G_beta,
    double rho    = 1.0,
    int num_threads = 1) {

  int n = D.n_rows;
  int p = G_beta.n_rows;
  int T = Y.n_slices;
  int J = par_grid.n_rows;

  bool G_is_id = (rho == 1.0) &&
    arma::approx_equal(G_beta, arma::eye<arma::mat>(p, p),
                       "absdiff", 1e-12);

  // -- Precompute K_vec and iK_vec once -------------
  std::vector<arma::mat> K_vec(J), iK_vec(J);
  std::vector<double>    c_vec(J);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(num_threads)
#endif
  for (int j = 0; j < J; j++) {
    K_vec[j] = arma::exp(-par_grid(j, 1) * D);
    if (!arma::inv_sympd(iK_vec[j], K_vec[j])) {
      arma::mat Kj_reg = K_vec[j];
      Kj_reg.diag() += 1e-10;
      arma::inv_sympd(iK_vec[j], Kj_reg);
    }
    c_vec[j] = (1.0 - par_grid(j, 0)) / par_grid(j, 0);
  }

  // -- Replicate prior J times --------------------------------------
  std::vector<Rcpp::List> priors_vec(J, prior);
  Rcpp::List FF_prec = Rcpp::wrap(priors_vec);
  Rcpp::List out_FF(T);

  // -- Time loop ----------------------------------------------------
  // State carried between iterations as C++ vectors (no R list overhead)
  std::vector<arma::mat> m_vec(J), C_bb_vec(J), C_bw_vec(J), C_ww_vec(J);
  std::vector<arma::mat> Psi_vec(J);
  std::vector<double>    nu_vec(J);

  for (int j = 0; j < J; j++) {
    Rcpp::List p_j = FF_prec[j];
    m_vec[j]    = Rcpp::as<arma::mat>(p_j["m"]);
    C_bb_vec[j] = Rcpp::as<arma::mat>(p_j["C_bb"]);
    C_bw_vec[j] = Rcpp::as<arma::mat>(p_j["C_bw"]);
    C_ww_vec[j] = Rcpp::as<arma::mat>(p_j["C_ww"]);
    Psi_vec[j]  = Rcpp::as<arma::mat>(p_j["Psi"]);
    nu_vec[j]   = Rcpp::as<double>(p_j["nu"]);
  }

  for (int t = 0; t < T; t++) {
    arma::mat Yt = Y.slice(t);
    arma::mat Xt = X.slice(t);

    std::vector<FFResult_v2> results;
    arma::mat density_evaluations;
    run_filters_and_density(Yt, Xt, G_beta, rho, G_is_id,
                            K_vec, c_vec,
                            m_vec, C_bb_vec, C_bw_vec, C_ww_vec,
                            nu_vec, Psi_vec,
                            J, n, num_threads,
                            results, density_evaluations);

    // Update state vectors in-place (avoids R list round-trip)
    for (int j = 0; j < J; j++) {
      m_vec[j]    = results[j].m;
      C_bb_vec[j] = results[j].C_bb;
      C_bw_vec[j] = results[j].C_bw;
      C_ww_vec[j] = results[j].C_ww;
      Psi_vec[j]  = results[j].Psi;
      nu_vec[j]   = results[j].nu;
    }

    // Pack into R list for out_FF (needed by backward sampler/R)
    Rcpp::List filtered_results(J);
    for (int j = 0; j < J; j++) {
      filtered_results[j] = Rcpp::List::create(
        Rcpp::Named("m")    = results[j].m,
        Rcpp::Named("C_bb") = results[j].C_bb,
        Rcpp::Named("C_bw") = results[j].C_bw,
        Rcpp::Named("C_ww") = results[j].C_ww,
        Rcpp::Named("nu")   = results[j].nu,
        Rcpp::Named("Psi")  = results[j].Psi);
    }
    out_FF[t] = Rcpp::List::create(
      Rcpp::Named("filtered_results")    = filtered_results,
      Rcpp::Named("density_evaluations") = density_evaluations);
  }
  return out_FF;
}


// ============================================================
// SECTION 9 - sample_index_v2
// ============================================================

// [[Rcpp::export]]
arma::uvec sample_index_v2(const int& size, const int& length,
                           const arma::vec& p) {
  arma::uvec seq = arma::linspace<arma::uvec>(0, size - 1, size);
  return Rcpp::RcppArmadillo::sample(seq, length, true, p);
}


// ============================================================
// SECTION 10 - Block Schur backward cache
// ============================================================

struct BackwardCache {
  arma::mat H;       // np x np
  arma::mat L_H;     // np x np  lower Chol of H
  arma::mat eta;     // np x q   = iC * m (block form)
  arma::mat iK_j;    // n  x n
  int p_dim, n_dim;
};

// Fills a pre-allocated BackwardCache using standard Armadillo operations.
// Called SERIALLY - no OMP inside, no BLAS over-subscription risk.
static void build_backward_cache(
    BackwardCache& bc,
    const arma::mat& m,
    const arma::mat& C_bb,
    const arma::mat& C_bw,
    const arma::mat& C_ww,
    const arma::mat& iK_j,
    const arma::mat& G_beta,
    double rho,
    bool G_is_identity) {

  int p = static_cast<int>(C_bb.n_rows);
  int n = static_cast<int>(C_ww.n_rows);
  bc.p_dim = p;  bc.n_dim = n;

  // iC_bb: p x p - arma::inv_sympd safe
  arma::mat iC_bb;
  if (!arma::inv_sympd(iC_bb, C_bb)) {
    arma::mat r = C_bb;
    r.diag() += 1e-10 * arma::trace(C_bb) / p;
    arma::inv_sympd(iC_bb, r);
  }

  // Schur complement S = C_ww - C_bw' * iC_bb * C_bw
  arma::mat T1 = iC_bb * C_bw;
  arma::mat S  = C_ww - C_bw.t() * T1;
  S = 0.5 * (S + S.t());

  // inv(S): n x n - arma::inv_sympd.
  arma::mat iS;
  if (!arma::inv_sympd(iS, S)) {
    S.diag() += 1e-10 * arma::trace(S) / n;
    arma::inv_sympd(iS, S);
  }

  // iC blocks
  arma::mat iC_bw  = -(T1 * iS);
  arma::mat iC_bb2 = iC_bb + (T1 * iS) * T1.t();

  // G'*iW_j*G blocks
  arma::mat GtiWG_bb, GtiWG_ww;
  if (G_is_identity) {
    GtiWG_bb = arma::eye<arma::mat>(p, p);
    GtiWG_ww = iK_j;
  } else {
    GtiWG_bb = G_beta.t() * G_beta;
    GtiWG_ww = (rho * rho) * iK_j;
  }

  // H^{-1} blocks
  arma::mat Hi_bb = iC_bb2 + GtiWG_bb;
  arma::mat Hi_bw = iC_bw;
  arma::mat Hi_ww = iS + GtiWG_ww;

  // iHi_bb: p x p - arma safe
  arma::mat iHi_bb;
  if (!arma::inv_sympd(iHi_bb, Hi_bb)) {
    Hi_bb.diag() += 1e-10;
    arma::inv_sympd(iHi_bb, Hi_bb);
  }
  arma::mat T2  = Hi_bw.t() * iHi_bb;
  arma::mat S_H = Hi_ww - T2 * Hi_bw;
  S_H = 0.5 * (S_H + S_H.t());

  // inv(S_H): n x n - arma
  arma::mat iS_H;
  if (!arma::inv_sympd(iS_H, S_H)) {
    S_H.diag() += 1e-10 * arma::trace(S_H) / n;
    arma::inv_sympd(iS_H, S_H);
  }

  arma::mat H_bw = -(iHi_bb * Hi_bw) * iS_H;
  arma::mat H_bb = iHi_bb + iHi_bb * Hi_bw * iS_H * T2;

  // Assemble H into pre-allocated np x np matrix
  bc.H.submat(0, 0, p-1, p-1)     = 0.5 * (H_bb + H_bb.t());
  bc.H.submat(0, p, p-1, p+n-1)   = H_bw;
  bc.H.submat(p, 0, p+n-1, p-1)   = H_bw.t();
  bc.H.submat(p, p, p+n-1, p+n-1) = 0.5 * (iS_H + iS_H.t());

  // Cholesky of H - safe_chol handles jitter fallback
  arma::mat H_sym = 0.5 * (bc.H + bc.H.t());
  bc.L_H = safe_chol(H_sym);

  // eta = iC * m in block form
  arma::mat m_b = m.rows(0, p-1);
  arma::mat m_w = m.rows(p, p+n-1);
  bc.eta.rows(0, p-1)   = iC_bb2 * m_b + iC_bw * m_w;
  bc.eta.rows(p, p+n-1) = iC_bw.t() * m_b + iS * m_w;

  bc.iK_j = iK_j;
}

// [[Rcpp::export]]
arma::cube weighted_backward_sample_v2(
    const arma::mat& G_beta,
    double rho,
    bool G_is_identity,
    const arma::mat& D,
    const Rcpp::List& FF_t,
    const arma::cube& ThetaSmp,
    const arma::cube& SigmaSmp,
    const arma::mat&  par_grid,
    const arma::vec&  weights,
    int num_threads = 1,
    SEXP iK_vec_sexp = R_NilValue) {

  int n  = D.n_rows;
  int p  = G_beta.n_rows;
  int np = p + n;
  int J  = par_grid.n_rows;
  int L  = ThetaSmp.n_slices;
  int q  = ThetaSmp.n_cols;

  arma::uvec model_idx = sample_index_v2(J, L, weights);
  arma::uvec uniqmod   = arma::unique(model_idx);
  int        JJ        = static_cast<int>(uniqmod.n_elem);

  std::unordered_map<arma::uword, arma::uword> mod_to_local;
  mod_to_local.reserve(JJ);
  for (int j = 0; j < JJ; j++)
    mod_to_local[uniqmod(j)] = static_cast<arma::uword>(j);

  // iK_vec: use precomputed if passed, else compute
  std::vector<arma::mat> iK_local(J);
  if (!Rf_isNull(iK_vec_sexp)) {
    Rcpp::List iKl(iK_vec_sexp);
    for (int j = 0; j < J; j++)
      iK_local[j] = Rcpp::as<arma::mat>(iKl[j]);
  } else {
    for (int j = 0; j < J; j++) {
      arma::mat K_j = arma::exp(-par_grid(j, 1) * D);
      arma::inv_sympd(iK_local[j], K_j);
    }
  }

  // -- Build backward cache per unique model ----------------------------
  std::vector<BackwardCache> cache(JJ);
  for (int jj = 0; jj < JJ; jj++) {
    cache[jj].H.zeros(np, np);    // pre-allocate serially
    cache[jj].L_H.zeros(np, np);
    cache[jj].eta.zeros(np, q);
    cache[jj].iK_j.zeros(n, n);
    cache[jj].p_dim = p;
    cache[jj].n_dim = n;
  }
  // Phase 1 (serial): Rcpp extraction - R API is not thread-safe
  std::vector<arma::mat> m_ex(JJ), Cbb_ex(JJ), Cbw_ex(JJ), Cww_ex(JJ);
  for (int jj = 0; jj < JJ; jj++) {
    arma::uword j_mod = uniqmod(jj);
    Rcpp::List FF_j   = FF_t(j_mod);
    m_ex[jj]   = Rcpp::as<arma::mat>(FF_j["m"]);
    Cbb_ex[jj] = Rcpp::as<arma::mat>(FF_j["C_bb"]);
    Cbw_ex[jj] = Rcpp::as<arma::mat>(FF_j["C_bw"]);
    Cww_ex[jj] = Rcpp::as<arma::mat>(FF_j["C_ww"]);
  }

  // Phase 2 (parallel): fill cache - noblas functions - no BLAS thread pool
  // Each thread calls inv_sympd_noblas / chol_lower_noblas: pure C++ loops,
  // zero global state, zero OpenBLAS interaction - no over-subscription.
  int cache_threads = std::min(num_threads, JJ);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(cache_threads)
#endif
  for (int jj = 0; jj < JJ; jj++) {
    arma::uword j_mod = uniqmod(jj);
    build_backward_cache(
      cache[jj],
           m_ex[jj], Cbb_ex[jj], Cbw_ex[jj], Cww_ex[jj],
                                                   iK_local[j_mod], G_beta, rho, G_is_identity);
  }

  // Draw L samples in parallel
  arma::cube Theta_ss(np, q, L);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(num_threads)
#endif
  for (int l = 0; l < L; l++) {
    arma::uword l_mod = model_idx(l);
    arma::uword lj    = mod_to_local.at(l_mod);
    const BackwardCache& bc = cache[lj];
    int p_b = bc.p_dim, n_b = bc.n_dim;

    arma::mat Theta_l = ThetaSmp.slice(l);
    arma::mat Sigma_l = SigmaSmp.slice(l);

    // iW_j * Theta_l in block form
    arma::mat iWT(np, q);
    iWT.rows(0, p_b-1)      = Theta_l.rows(0, p_b-1);
    iWT.rows(p_b, p_b+n_b-1) = bc.iK_j * Theta_l.rows(p_b, p_b+n_b-1);

    arma::mat GiWT(np, q);
    if (G_is_identity) {
      GiWT = iWT;
    } else {
      GiWT.rows(0, p_b-1)      = G_beta.t() * iWT.rows(0, p_b-1);
      GiWT.rows(p_b, p_b+n_b-1) = rho * iWT.rows(p_b, p_b+n_b-1);
    }

    arma::mat h    = bc.H * (bc.eta + GiWT);
    arma::mat L_Sg = safe_chol(Sigma_l);
    arma::mat Z    = std_normal_mat(np, q);
    Theta_ss.slice(l) = h + bc.L_H * Z * L_Sg.t();
  }
  return Theta_ss;
}


// ============================================================
// SECTION 11 - Multi-step backward sampling
// ============================================================

// [[Rcpp::export]]
Rcpp::List weighted_backward_sample_T_v2(
    const arma::mat&  G_beta,
    double rho,
    const arma::mat&  D,
    const Rcpp::List& ForwFilt,
    const int L,
    const arma::mat&  par_grid,
    const arma::vec&  weights,
    int num_threads = 1) {

  int J    = par_grid.n_rows;
  int Tmax = static_cast<int>(ForwFilt.size()) - 1;
  bool G_is_id = (rho == 1.0) &&
    arma::approx_equal(
      G_beta, arma::eye<arma::mat>(G_beta.n_rows, G_beta.n_rows),
      "absdiff", 1e-12);

  // Precompute iK_vec once for the entire backward pass
  std::vector<arma::mat> iK_vec(J);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(num_threads)
#endif
  for (int j = 0; j < J; j++) {
    arma::mat K_j = arma::exp(-par_grid(j, 1) * D);
    arma::inv_sympd(iK_vec[j], K_j);
  }
  Rcpp::List iK_list(J);
  for (int j = 0; j < J; j++) iK_list[j] = iK_vec[j];
  SEXP iK_sexp = Rcpp::wrap(iK_list);

  Rcpp::List out_BS(Tmax + 1);

  // -- Terminal step: MNIW draw --------------------
  Rcpp::List out_FFtmax = ForwFilt(Tmax);
  Rcpp::List FF_T       = out_FFtmax["filtered_results"];

  arma::uvec last_idx = sample_index_v2(J, L, weights);
  arma::uvec uniqmod  = arma::unique(last_idx);
  int JJ = static_cast<int>(uniqmod.n_elem);

  arma::uvec LLi(JJ);
  for (int j = 0; j < JJ; j++)
    LLi(j) = arma::accu(last_idx == uniqmod(j));

  Rcpp::List FF_j0 = FF_T(uniqmod(0));
  arma::mat m0 = Rcpp::as<arma::mat>(FF_j0["m"]);
  int np = m0.n_rows, q = m0.n_cols, p = G_beta.n_rows;
  arma::cube ThetaSmpT(np, q, L);
  arma::cube SigmaSmpT(q,  q, L);

  int slice_start = 0;
  for (int j = 0; j < JJ; j++) {
    arma::uword j_mod  = uniqmod(j);
    Rcpp::List FF_j    = FF_T(j_mod);
    arma::mat m_j      = Rcpp::as<arma::mat>(FF_j["m"]);
    arma::mat C_bb_j   = Rcpp::as<arma::mat>(FF_j["C_bb"]);
    arma::mat C_bw_j   = Rcpp::as<arma::mat>(FF_j["C_bw"]);
    arma::mat C_ww_j   = Rcpp::as<arma::mat>(FF_j["C_ww"]);
    double    nu_j     = Rcpp::as<double>(FF_j["nu"]);
    arma::mat Psi_j    = Rcpp::as<arma::mat>(FF_j["Psi"]);

    // Build C_full once (terminal step, only JJ times total)
    int n_j = C_ww_j.n_rows;
    arma::mat C_full(np, np);
    C_full.submat(0, 0, p-1, p-1)         = C_bb_j;
    C_full.submat(0, p, p-1, p+n_j-1)     = C_bw_j;
    C_full.submat(p, 0, p+n_j-1, p-1)     = C_bw_j.t();
    C_full.submat(p, p, p+n_j-1, p+n_j-1) = C_ww_j;
    arma::mat L_C = safe_chol(C_full);

    int LL = static_cast<int>(LLi(j));
    arma::cube Th_j, Si_j;
    rMNIW_batch_Lc(LL, m_j, L_C, nu_j, Psi_j, Th_j, Si_j, num_threads);

    ThetaSmpT.slices(slice_start, slice_start + LL - 1) = Th_j;
    SigmaSmpT.slices(slice_start, slice_start + LL - 1) = Si_j;
    slice_start += LL;
  }

  arma::uvec shuf         = arma::randperm(L);
  arma::cube ThetaSmpTmax = ThetaSmpT.slices(shuf);
  arma::cube SigmaSmpTmax = SigmaSmpT.slices(shuf);

  out_BS(Tmax) = weighted_backward_sample_v2(
    G_beta, rho, G_is_id, D,
    FF_T, ThetaSmpTmax, SigmaSmpTmax,
    par_grid, weights, num_threads, iK_sexp);

  arma::cube SigmaCarry = SigmaSmpTmax;

  for (int t = Tmax - 1; t >= 0; t--) {
    Rcpp::List ForwFilt_t = ForwFilt(t);
    Rcpp::List out_FF_t   = ForwFilt_t["filtered_results"];
    arma::cube ThetaSmpBS = out_BS(t + 1);
    out_BS(t) = weighted_backward_sample_v2(
      G_beta, rho, G_is_id, D,
      out_FF_t, ThetaSmpBS, SigmaCarry,
      par_grid, weights, num_threads, iK_sexp);
  }
  return out_BS;
}


// ============================================================
// SECTION 12 - Weight computation
// ============================================================

// Internal pure-C++ implementation: works on raw double arrays only.
// Safe to call from OMP threads - no Rcpp, no R heap access.
static std::vector<double> weights_proj_impl(
    const arma::mat& scores,   // (t-1) x J  pre-filled, log-scale exped
    double lr, int max_iter) {
  int n_obs = static_cast<int>(scores.n_rows);
  int J     = static_cast<int>(scores.n_cols);
  std::vector<double> w(J, 1.0 / J);
  for (int iter = 0; iter < max_iter; ++iter) {
    std::vector<double> grad(J, 0.0);
    for (int i = 0; i < n_obs; ++i) {
      double dot = 0.0;
      for (int j = 0; j < J; ++j) dot += scores(i, j) * w[j];
      if (dot <= 0.0) continue;
      for (int j = 0; j < J; ++j) grad[j] += scores(i, j) / dot;
    }
    for (int j = 0; j < J; ++j) { grad[j] /= n_obs; w[j] += lr * grad[j]; }
    int nw = J;
    std::vector<double> u(w);
    std::sort(u.begin(), u.end(), std::greater<double>());
    double cumsum = 0.0, rho_s = 0.0;
    for (int j = 0; j < nw; ++j) {
      cumsum += u[j];
      if (u[j] > (cumsum - 1.0) / (j + 1)) rho_s = j + 1;
    }
    double theta = (std::accumulate(u.begin(),
                                    u.begin() + static_cast<int>(rho_s),
                                    0.0) - 1.0) / rho_s;
    for (int j = 0; j < nw; ++j) w[j] = std::max(w[j] - theta, 0.0);
  }
  return w;
}

// Exported wrapper: keeps original R-facing signature unchanged.
// [[Rcpp::export]]
Rcpp::NumericVector optimize_weights_proj_v2(
    Rcpp::NumericMatrix scores, double lr = 0.05, int max_iter = 500) {
  int nr = scores.nrow(), nc = scores.ncol();
  arma::mat sc(nr, nc);
  for (int j = 0; j < nc; ++j)
    for (int i = 0; i < nr; ++i) sc(i, j) = scores(i, j);
  std::vector<double> w = weights_proj_impl(sc, lr, max_iter);
  return Rcpp::NumericVector(w.begin(), w.end());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_Wt_cpp_v2(
    Rcpp::List density_list, int n, int t,
    double lr = 0.05, int max_iter = 500, int n_threads = 1) {

  // -- SERIAL pre-extraction: ALL Rcpp access happens here, before OMP ---
  // Rcpp objects must never be created or accessed inside an OMP region -
  // R's heap is not thread-safe.  We extract everything to arma::mat now.
  int T1 = t - 1;   // number of density slices
  int J  = Rcpp::as<arma::mat>(density_list[0]).n_cols;

  // dens[a](i,j) = exp(log-density at location i, model j, time a+1)
  std::vector<arma::mat> dens(T1);
  for (int a = 0; a < T1; ++a) {
    arma::mat raw = Rcpp::as<arma::mat>(density_list[a]);  // n x J log-dens
    dens[a] = arma::exp(raw);                              // exponentiate
  }

  // Pre-allocate output as arma::mat (no R heap inside OMP)
  arma::mat Wi_arma(n, J);

  // -- PARALLEL loop over locations - pure arma/C++, zero Rcpp ----------
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n; ++i) {
    // Build (T1 x J) score matrix for location i
    arma::mat epd_i(T1, J);
    for (int a = 0; a < T1; ++a) {
      for (int j = 0; j < J; ++j)
        epd_i(a, j) = dens[a](i, j);
    }
    // Optimise weights (pure C++, no Rcpp)
    std::vector<double> w = weights_proj_impl(epd_i, lr, max_iter);
    for (int j = 0; j < J; ++j) Wi_arma(i, j) = w[j];
  }

  // -- SERIAL post-copy: arma - Rcpp return value -----------------------
  Rcpp::NumericMatrix Wi(n, J);
  for (int j = 0; j < J; ++j)
    for (int i = 0; i < n; ++i)
      Wi(i, j) = Wi_arma(i, j);
  return Wi;
}


// ============================================================
// SECTION 13 - Posterior predictive
// ============================================================

// [[Rcpp::export]]
Rcpp::List predict_ffbs_v2(
    const arma::mat& G_beta,
    double rho,
    bool G_is_identity,
    const arma::mat& X_new,
    const arma::mat& K_j,
    double c_j,
    const arma::mat& a,
    const arma::mat& C_bb, const arma::mat& C_bw, const arma::mat& C_ww,
    double nu, const arma::mat& Psi,
    int L, int num_threads = 1) {

  int np = a.n_rows, q = a.n_cols, p = G_beta.n_rows, n = np - p;

  arma::mat a_new(np, q);
  if (G_is_identity) {
    a_new = a;
  } else {
    a_new.rows(0, p-1)  = G_beta * a.rows(0, p-1);
    a_new.rows(p, np-1) = rho * a.rows(p, np-1);
  }

  arma::mat R_bb, R_bw, R_ww;
  if (G_is_identity) {
    R_bb = C_bb; R_bb.diag() += 1.0;
    R_bw = C_bw;
    R_ww = C_ww + K_j;
  } else {
    R_bb = G_beta * C_bb * G_beta.t(); R_bb.diag() += 1.0;
    R_bw = rho * (G_beta * C_bw);
    R_ww = (rho*rho)*C_ww + K_j;
  }

  arma::mat f_pred = X_new * a_new.rows(0, p-1) + a_new.rows(p, np-1);
  arma::mat A_top  = R_bb * X_new.t() + R_bw;
  arma::mat Q_pred = R_ww + X_new * A_top;
  Q_pred.diag() += c_j;
  Q_pred = 0.5*(Q_pred + Q_pred.t());

  arma::mat L_Q;
  if (!arma::chol(L_Q, Q_pred, "lower")) L_Q = safe_chol(Q_pred);

  arma::cube Y_pred(n, q, L);
  spFFBS_par_for(L, num_threads, [&](int l) {
    arma::mat Sig = rIW(nu, Psi);
    arma::mat L_S = safe_chol(Sig);
    arma::mat Z   = std_normal_mat(n, q);
    Y_pred.slice(l) = f_pred + L_Q * Z * L_S.t();
  });

  return Rcpp::List::create(
    Rcpp::Named("Y")    = Y_pred,
    Rcpp::Named("a")    = a_new,
    Rcpp::Named("R_bb") = R_bb,
    Rcpp::Named("R_bw") = R_bw,
    Rcpp::Named("R_ww") = R_ww);
}


// ============================================================
// SECTION 14 - Temporal forecast
// ============================================================

// [[Rcpp::export]]
Rcpp::List temporal_forecast_v2(
    const arma::mat&  G_beta,
    double rho,
    const arma::mat&  D,
    const arma::mat&  par_grid,
    const Rcpp::List& ForwFilt,
    const arma::cube& X_all,
    const arma::mat&  weights,
    int horiz, int L, int num_threads = 1) {

  int p    = G_beta.n_rows;
  int J    = par_grid.n_rows;
  int tmax = weights.n_cols + 1;
  int th   = tmax + horiz;

  bool G_is_id = (rho == 1.0) &&
    arma::approx_equal(G_beta, arma::eye<arma::mat>(p, p),
                       "absdiff", 1e-12);

  // Precompute K_vec and c_vec
  std::vector<arma::mat> K_vec(J);
  std::vector<double>    c_vec(J);
  for (int j = 0; j < J; j++) {
    K_vec[j] = arma::exp(-par_grid(j, 1) * D);
    c_vec[j] = (1.0 - par_grid(j, 0)) / par_grid(j, 0);
  }

  Rcpp::List pred(th), predictions(th);

  for (int t = 0; t < th; t++) {
    arma::uword w_idx = static_cast<arma::uword>(std::min(t, tmax - 1));
    arma::vec weights_t;
    if (t < 1) {
      weights_t = arma::vec(J, arma::fill::ones) / static_cast<double>(J);
    } else {
      weights_t = arma::vec(weights.col(w_idx - 1));
    }
    arma::uvec  midx = sample_index_v2(J, 1, weights_t);
    arma::uword j    = midx(0);
    arma::mat   X_t  = X_all.slice(t);

    if (t < tmax) {
      Rcpp::List FF_t = ForwFilt(t);
      Rcpp::List FF_J = FF_t["filtered_results"];
      Rcpp::List FF_j = FF_J(j);
      Rcpp::List p_t  = predict_ffbs_v2(
        G_beta, rho, G_is_id, X_t, K_vec[j], c_vec[j],
                                                  Rcpp::as<arma::mat>(FF_j["m"]),
                                                  Rcpp::as<arma::mat>(FF_j["C_bb"]),
                                                  Rcpp::as<arma::mat>(FF_j["C_bw"]),
                                                  Rcpp::as<arma::mat>(FF_j["C_ww"]),
                                                  Rcpp::as<double>(FF_j["nu"]),
                                                  Rcpp::as<arma::mat>(FF_j["Psi"]),
                                                  L, num_threads);
      pred[t]        = p_t;
      predictions[t] = p_t["Y"];
    } else {
      Rcpp::List prev_p  = pred[t - 1];
      Rcpp::List FF_last = ForwFilt(tmax - 1);
      Rcpp::List FF_J_l  = FF_last["filtered_results"];
      Rcpp::List FF_j_l  = FF_J_l(j);
      Rcpp::List p_t     = predict_ffbs_v2(
        G_beta, rho, G_is_id, X_t, K_vec[j], c_vec[j],
                                                  Rcpp::as<arma::mat>(prev_p["a"]),
                                                  Rcpp::as<arma::mat>(prev_p["R_bb"]),
                                                  Rcpp::as<arma::mat>(prev_p["R_bw"]),
                                                  Rcpp::as<arma::mat>(prev_p["R_ww"]),
                                                  Rcpp::as<double>(FF_j_l["nu"]),
                                                  Rcpp::as<arma::mat>(FF_j_l["Psi"]),
                                                  L, num_threads);
      pred[t]        = p_t;
      predictions[t] = p_t["Y"];
    }
  }
  return predictions;
}


// ============================================================
// SECTION 15 - Spatial interpolation
// ============================================================

// [[Rcpp::export]]
Rcpp::List spatial_interpolation_v2(
    const arma::mat&  G_beta,
    double rho,
    const arma::mat&  Xu,
    const arma::mat&  D_s,
    const arma::mat&  D_u,
    const arma::mat&  D_us,
    const arma::mat&  par_grid,
    const Rcpp::List& ForwFilt,
    const arma::mat&  weights,
    int t, int L, int num_threads = 1) {

  int u  = D_u.n_rows;
  int p  = G_beta.n_rows;
  int J  = par_grid.n_rows;
  int T_obs = static_cast<int>(ForwFilt.size());  // observed time horizon

  // -- Support t > T_obs (out-of-sample spatial prediction)  -------------
  // For t <= T_obs: use filtered posterior at time t.
  // For t  > T_obs: propagate the state k = t - T_obs steps beyond T_obs
  //   using the random-walk state equation (G = I):
  //     m_{T+k} = m_T                    (mean unchanged under RW)
  //     C_bb,{T+k} ≈ C_bb,T              (beta covariance, W_B negligible)
  //     C_bw,{T+k}  = C_bw,T             (cross-block unchanged)
  //     C_ww,{T+k}  = C_ww,T + k * K_j  (spatial RW accumulates k innovations)
  // The Wtilde residual additionally carries the factor (T_obs + k) instead of t.
  int t_ref    = std::min(t, T_obs);           // index into ForwFilt (1-based)
  int k_fcast  = std::max(0, t - T_obs);       // out-of-sample steps (0 if in-sample)

  // Weight vector: use last column if t > number of weight columns
  arma::vec w0 = arma::vec(J, arma::fill::ones) / J;
  arma::mat wmat = arma::join_horiz(w0, weights);  // J x (ncol(weights)+1)
  int w_col = std::min(t - 1, static_cast<int>(wmat.n_cols) - 1);
  arma::vec weights_t = wmat.col(w_col);

  arma::uvec model_idx = sample_index_v2(J, L, weights_t);
  arma::uvec uniqmod   = arma::unique(model_idx);
  int JJ = static_cast<int>(uniqmod.n_elem);

  std::unordered_map<arma::uword, arma::uword> mod_to_local;
  mod_to_local.reserve(JJ);
  for (int j = 0; j < JJ; j++)
    mod_to_local[uniqmod(j)] = static_cast<arma::uword>(j);

  struct SpatialModel {
    arma::mat mu; arma::mat L_E; double nu; arma::mat Psi;
    int rows, cols;
  };
  std::vector<SpatialModel> models(JJ);

  for (int jj = 0; jj < JJ; jj++) {
    arma::uword j_mod = uniqmod(jj);
    double tau_j = par_grid(j_mod, 0);
    double phi_j = par_grid(j_mod, 1);

    arma::mat Vtilde = ((1.0 - tau_j) / tau_j) * arma::eye<arma::mat>(u, u);
    arma::mat Kss    = arma::exp(-phi_j * D_s);
    arma::mat Kuu    = arma::exp(-phi_j * D_u);
    arma::mat Kus    = arma::exp(-phi_j * D_us);
    arma::mat iKss;
    arma::inv_sympd(iKss, Kss);
    arma::mat Mtilde = Kus * iKss;

    // W_tilde scales linearly with t
    arma::mat Wtilde = static_cast<double>(t) * (Kuu - Mtilde * Kus.t());
    Wtilde = 0.5*(Wtilde + Wtilde.t());

    // Extract filtered state at t_ref (last observed time if out-of-sample)
    Rcpp::List FF_t  = ForwFilt(t_ref - 1);
    Rcpp::List FF_J  = FF_t["filtered_results"];
    Rcpp::List FF_j  = FF_J(j_mod);
    arma::mat C_bb_j = Rcpp::as<arma::mat>(FF_j["C_bb"]);
    arma::mat C_bw_j = Rcpp::as<arma::mat>(FF_j["C_bw"]);
    arma::mat C_ww_j = Rcpp::as<arma::mat>(FF_j["C_ww"]);
    arma::mat m_j    = Rcpp::as<arma::mat>(FF_j["m"]);
    arma::mat Psi_j  = Rcpp::as<arma::mat>(FF_j["Psi"]);
    double    nu_j   = Rcpp::as<double>(FF_j["nu"]);
    int       n_j    = C_ww_j.n_rows;
    int       np_j   = p + n_j;

    // Propagate state k_fcast steps for out-of-sample t
    if (k_fcast > 0) {
      C_ww_j += static_cast<double>(k_fcast) * Kss;

    }

    arma::mat C_full(np_j, np_j);
    C_full.submat(0, 0, p-1, p-1)         = C_bb_j;
    C_full.submat(0, p, p-1, p+n_j-1)     = C_bw_j;
    C_full.submat(p, 0, p+n_j-1, p-1)     = C_bw_j.t();
    C_full.submat(p, p, p+n_j-1, p+n_j-1) = C_ww_j;

    arma::mat zero_up = arma::zeros<arma::mat>(u, p);
    arma::mat chi_t   = arma::join_horiz(
      arma::join_vert(Xu, zero_up),
      arma::join_vert(Mtilde, Mtilde));
    arma::mat N_t     = arma::join_horiz(
      arma::join_vert(Vtilde + Wtilde, Wtilde),
      arma::join_vert(Wtilde, Wtilde));

    arma::mat mu_t = chi_t * m_j;
    arma::mat E_t  = chi_t * C_full * chi_t.t() + N_t;
    E_t = 0.5*(E_t + E_t.t());
    arma::mat L_E = safe_chol(E_t);  // precomputed once

    models[jj] = {mu_t, L_E, nu_j, Psi_j,
                  static_cast<int>(mu_t.n_rows),
                  static_cast<int>(mu_t.n_cols)};
  }

  std::vector<arma::mat> preds(L);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(num_threads)
#endif
  for (int l = 0; l < L; l++) {
    arma::uword l_mod = model_idx(l);
    arma::uword lj    = mod_to_local.at(l_mod);
    const SpatialModel& sm = models[lj];
    arma::mat Sig = rIW(sm.nu, sm.Psi);
    arma::mat L_S = safe_chol(Sig);
    arma::mat Z   = std_normal_mat(sm.rows, sm.cols);
    preds[l] = sm.mu + sm.L_E * Z * L_S.t();
  }

  Rcpp::List predictions(L);
  for (int l = 0; l < L; l++) predictions[l] = preds[l];
  return predictions;
}
