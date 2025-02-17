#include <RcppArmadillo.h>
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// PARALLEL FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// FORWARD FILTERING
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward Filtering - SingleStep
//[[Rcpp::export(name = "forward_filter")]]
List forward_filter(const arma::mat& Y, const arma::mat& G, const arma::mat& P, const arma::mat& V, const arma::mat& W, const arma::mat& m, const arma::mat& C, const double& nu, const arma::mat& Psi) {

  // Filtered prior
  arma::mat a = G * m;
  arma::mat R = G * C * trans(G) + W;

  arma::mat P_trans = trans(P);
  // 1-Step ahead predictive
  arma::mat f = P * a;
  arma::mat Q = P * R * P_trans + V;

  // Linear Systems Avoiding Inverse
  arma::mat z = solve(V, P); // V^-1P
  arma::mat x = solve(V, Y); // V^-1Y
  arma::mat y = solve(Q, (Y-f));

  // Inverse Matrices
  arma::mat R_inv = inv(R);

  // Filtered posterior
  arma::mat C_inv = R_inv + (P_trans * z);
  arma::mat C_new = inv(C_inv);
  arma::mat m_new = C_new * (R_inv * a + (P_trans * x));;

  // Linear System Avoiding Inverse
  // Non-dynamic parameters posterior
  double nu_new  = nu  + (Y.n_rows/2);
  arma::mat Psi_new = Psi + (0.5) * ( trans(Y - f) * y );

  // Return FF parameters as an R list
  return List::create(Named("a") = a,
                      Named("R") = R,
                      Named("f") = f,
                      Named("Q") = Q,
                      Named("m") = m_new,
                      Named("C") = C_new,
                      Named("nu") = nu_new,
                      Named("Psi") = Psi_new,
                      Named("iC") = C_inv);

}


// [[Rcpp::export]]
List parallel_forward_filter(const arma::mat& Y, const arma::mat& G, const arma::mat& P, const arma::mat& D,
                             const arma::mat& par_grid, List FF_prec, int J, int p, int n,
                             int num_threads = 1) {

  // Create cointainer
  List results(J);

  // Set the number of threads for OpenMP
  omp_set_num_threads(num_threads);

#pragma omp parallel for schedule(dynamic)
  for (int j = 0; j < J; j++) {
    try {
      // Extract parameters
      double phi_j = par_grid(j, 1);
      arma::mat K_j = exp(-phi_j * D);

      // Fix dimension mismatch in W_j
      arma::uword n_uword = static_cast<arma::uword>(n);

      // Rcpp::Rcout << "K_j dimensions: " << K_j.n_rows << " x " << K_j.n_cols << std::endl;
      // Rcpp::Rcout << "D dimensions: " << D.n_rows << " x " << D.n_cols << std::endl;
      // Rcpp::Rcout << "phi_j dimensions: " << phi_j << std::endl;
      // Rcpp::Rcout << "n dimensions: " << n_uword << std::endl;

      if (K_j.n_rows != n_uword || K_j.n_cols != n_uword) {
        throw std::runtime_error("Dimension mismatch: K_j must be (n, n)");
      }

      arma::mat W_j = join_cols(
        join_rows(eye(p, p), zeros(p, n)),
        join_rows(zeros(n, p), K_j)
      );

      double tau_j = par_grid(j, 0);
      arma::mat V_j = ((1 - tau_j) / tau_j) * eye(n, n);

      // Extract previous results safely
      List res_jprec;

#pragma omp critical
{
  res_jprec = as<List>(FF_prec[j]);
}

arma::mat m = as<arma::mat>(res_jprec["m"]);
arma::mat C = as<arma::mat>(res_jprec["C"]);
double nu = as<double>(res_jprec["nu"]);
arma::mat Psi = as<arma::mat>(res_jprec["Psi"]);

// Forward filtering
List out_j_t = forward_filter(Y, G, P, V_j, W_j, m, C, nu, Psi);

// Store results safely
#pragma omp critical
{
  results[j] = out_j_t;
  // results[j] = List::create(Named("out_j_t") = out_j_t);
}

    } catch (std::exception &ex) {
      Rcpp::Rcerr << "Exception in thread " << j << ": " << ex.what() << std::endl;
    }
  }

  return results;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// MATRIX T EVALUATIONS
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//' Multivariate Gamma function
//'
//' @export
// [[Rcpp::export]]
double lmvgamma(double x, int p) {
  double result = 0.25 * p * (p - 1) * log(M_PI);
  for (int j = 1; j <= p; ++j) {
    result += lgamma(x + 0.5 * (1 - j));
  }
  return result;
}


// Log-density of the Matrix T-distribution
// [[Rcpp::export]]
double log_matrix_t_density(const arma::mat& X, const arma::mat& M, const arma::mat& U, const arma::mat& V, double nu) {
  int n = X.n_rows;
  int p = X.n_cols;

  // Scale U and V by nu (to match mniw parameterization)
  arma::mat U_scaled = U; // / nu;
  arma::mat V_scaled = V; // / nu;

  // Compute Cholesky decompositions of U_scaled and V_scaled
  arma::mat L_U = chol(U_scaled, "lower");
  arma::mat L_V = chol(V_scaled, "lower");

  // Compute log-determinants of U_scaled and V_scaled
  double log_det_U = 2 * sum(log(L_U.diag()));
  double log_det_V = 2 * sum(log(L_V.diag()));

  // Compute the difference matrix
  arma::mat diff = X - M;

  // // Scale the difference matrix by the Cholesky factors
  // arma::mat iU = arma::inv(U_scaled);
  // arma::mat iV = arma::inv(V_scaled);
  // arma::mat Q = iU * diff * iV * diff.t();
  // double log_det_term = log_det(eye(n, n) + Q).real();

  // Scale the difference matrix by the Cholesky factors
  arma::mat scaled_diff = solve(trimatl(L_U), diff); // Solve L_U * Y = diff
  scaled_diff = solve(trimatl(L_V), scaled_diff.t()).t(); // Solve L_V * Z = Y^T, then transpose back
  // Compute the quadratic form Q = Z^T Z
  arma::mat Q = scaled_diff.t() * scaled_diff;
  // Compute the determinant term: log|I_n + Q|
  double log_det_term = log_det(eye(p, p) + Q).real();

  // Compute the log-density
  double log_density = lmvgamma(0.5 * (nu + n + p - 1), p) -
    lmvgamma(0.5 * (nu + p - 1), p) -
    (0.5 * n * p * log(M_PI)) -
    (0.5 * p * log_det_U) -
    (0.5 * n * log_det_V) -
    (0.5 * (nu + n + p - 1) * log_det_term);

  return log_density;
}

// Parallel Matrix T evaluations
// [[Rcpp::export]]
arma::mat parallel_matrix_t_evaluations(const arma::cube& Y, const Rcpp::List& out_J,
                                        int t, int J, int n,
                                        int num_threads = 1) {

  // Set the number of OpenMP threads
  omp_set_num_threads(num_threads);

  // Initialize output matrix (n x J)
  arma::mat pd_t(n, J, fill::zeros);

  // Parallel for loop over J
#pragma omp parallel for
  for (int j = 0; j < J; j++) {
    try {
      // Ensure valid indices
      if (j >= out_J.size()) {
#pragma omp critical
        Rcpp::Rcerr << "Index out of bounds for j = " << j << std::endl;
        continue;
      }

      Rcpp::List res_j = out_J[j];   // Get out_J[[j]]
      if (t >= res_j.size()) {
#pragma omp critical
        Rcpp::Rcerr << "Index out of bounds for t = " << t << ", J = " << j << std::endl;
        continue;
      }
      if (t - 1 >= res_j.size()) {
#pragma omp critical
        Rcpp::Rcerr << "Index out of bounds for t-1 = " << t-1 << ", J = " << j << std::endl;
        continue;
      }

      // Extract elements
      Rcpp::List res_jprec = res_j[t - 1];
      res_j = res_j[t];

      arma::mat f = res_j["f"];
      arma::mat Q = res_j["Q"];
      arma::mat Psi = res_jprec["Psi"];
      double nu = res_jprec["nu"];

      // Compute pd_t[g, j] in parallel
      for (int g = 0; g < n; g++) {
        if (static_cast<arma::uword>(g) >= Y.n_rows) continue;
        if (static_cast<arma::uword>(g) >= Q.n_rows || static_cast<arma::uword>(g) >= Q.n_cols) continue;

        // Extract Y[g, , t]
        arma::mat Ygt = Y.slice(t).row(g);

        // Rcpp::Rcout << "Ygt dimensions: " << Ygt << std::endl;
        // Rcpp::Rcout << "f dimensions: " << f.row(g) << std::endl;

        // Compute density using the log_matrix_t_density function
        double density = log_matrix_t_density(Ygt, f.row(g),
                                              Q.submat(g, g, g, g), Psi, nu);

#pragma omp critical
{
  pd_t(g, j) = density;
}
      }

    } catch (std::exception &ex) {
#pragma omp critical
      Rcpp::Rcerr << "Exception in thread " << j << ": " << ex.what() << std::endl;
    }
  }

  return pd_t;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// CONVEX OPTIMIZATION
///////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////
// UNIFIED FUNCTION
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Unified Parallel Function
// [[Rcpp::export]]
List unified_parallel_function(const arma::mat& Y, const arma::mat& G, const arma::mat& P, const arma::mat& D,
                               const arma::mat& par_grid, List FF_prec, int J, int p, int n,
                               int num_threads = 1) {

  // Set the number of threads for OpenMP
  omp_set_num_threads(num_threads);

  // Thread-safe containers for results
  std::vector<List> filtered_results(J);
  arma::mat density_evaluations(n, J, fill::zeros);

#pragma omp parallel for schedule(dynamic)
  for (int j = 0; j < J; j++) {
    try {
      // Extract parameters
      double phi_j = par_grid(j, 1);
      arma::mat K_j = exp(-phi_j * D);

      arma::uword n_uword = static_cast<arma::uword>(n);

      if (K_j.n_rows != n_uword || K_j.n_cols != n_uword) {
        throw std::runtime_error("Dimension mismatch: K_j must be (n, n)");
      }

      arma::mat W_j = join_cols(
        join_rows(eye(p, p), zeros(p, n)),
        join_rows(zeros(n, p), K_j)
      );

      double tau_j = par_grid(j, 0);
      arma::mat V_j = ((1 - tau_j) / tau_j) * eye(n, n);

      // Extract previous results safely
      List res_jprec;

#pragma omp critical
{
  res_jprec = as<List>(FF_prec[j]);
}

arma::mat m = as<arma::mat>(res_jprec["m"]);
arma::mat C = as<arma::mat>(res_jprec["C"]);
double nu = as<double>(res_jprec["nu"]);
arma::mat Psi = as<arma::mat>(res_jprec["Psi"]);

// Forward filtering
List out_j_t = forward_filter(Y, G, P, V_j, W_j, m, C, nu, Psi);

// Store filtered results
filtered_results[j] = out_j_t;

// Extract elements for density evaluation
arma::mat f = as<arma::mat>(out_j_t["f"]);
arma::mat Q = as<arma::mat>(out_j_t["Q"]);
arma::mat Psi_new = as<arma::mat>(out_j_t["Psi"]);
double nu_new = as<double>(out_j_t["nu"]);

// Compute density evaluations
for (int g = 0; g < n; g++) {
  if (static_cast<arma::uword>(g) >= Y.n_rows) continue;
  if (static_cast<arma::uword>(g) >= Q.n_rows || static_cast<arma::uword>(g) >= Q.n_cols) continue;

  arma::mat Ygt = Y.row(g);
  double density = log_matrix_t_density(Ygt, f.row(g),
                                        Q.submat(g, g, g, g), Psi_new, nu_new);

  density_evaluations(g, j) = density;
}

    } catch (std::exception &ex) {
#pragma omp critical
      Rcpp::Rcerr << "Exception in thread " << j << ": " << ex.what() << std::endl;
    }
  }

  // Return both filtered results and density evaluations
  return List::create(Named("filtered_results") = filtered_results,
                      Named("density_evaluations") = density_evaluations);
}


//' Parallel BPS-Forward Filtering - Multiple Step
//'
//' @export
// [[Rcpp::export]]
List spFF(const arma::cube& Y, const arma::mat& G, const arma::mat& P, const arma::mat& D,
          const arma::mat& par_grid, List const& prior,
          int num_threads = 1) {

  // Gather information from inputs
  int n = D.n_rows;
  int p = G.n_rows - n;
  int T = Y.n_slices;
  int J = par_grid.n_rows;

  // build containers
  List out_FF(T);

  // FirstStep - Prior Information
  // mat Y1 = Y.slice(0);
  // mat m0 = as<mat>(prior["m"]);
  // mat C0 = as<mat>(prior["C"]);
  // double nu0 = as<double>(prior["nu"]);
  // mat Psi0 = as<mat>(prior["Psi"]);
  // List out_prior = forward_filter(Y1, G, P, V, W, m0, C0, nu0, Psi0);
  arma::mat Y1 = Y.slice(0);;
  std::vector<List> priors = std::vector<List>(J, prior);
  List prior_list = Rcpp::wrap(priors);
  List out_prior = parallel_forward_filter(Y1, G, P, D, par_grid, prior_list, J, p, n, num_threads);
  out_FF(0) = List::create(Named("filtered_results") = out_prior);

  // OtherSteps - For Loop
  for (int t = 1; t < T; t++) {

    // Extract Old Information
    // mat Yt = Y.slice(t);
    // List out_old = out_FF(t-1);
    // mat m_old = as<mat>(out_old["m"]);
    // mat C_old = as<mat>(out_old["C"]);
    // double nu_old = as<double>(out_old["nu"]);
    // mat Psi_old = as<mat>(out_old["Psi"]);
    arma::mat Yt = Y.slice(t);
    List out_prec = out_FF(t-1);
    List out_old = out_prec["filtered_results"];

    // Compute New Information
    List out_new = unified_parallel_function(Yt, G, P, D, par_grid, out_old, J, p, n, num_threads);
    out_FF(t) = out_new;

  }

  // Return List of Results
  return out_FF;

}

