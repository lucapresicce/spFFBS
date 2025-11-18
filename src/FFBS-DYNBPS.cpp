#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// UTILITY FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//' Function to sample integers (index)
//'
//' @param size [integer] dimension of the set to sample
//' @param length [integer] number of elements to sample
//' @param p [vector] sampling probabilities
//'
//' @return [vector] sample of integers
//'
// [[Rcpp::export]]
arma::uvec sample_index(const int& size, const int& length, const arma::vec& p){
 arma::uvec sequence = arma::linspace<arma::uvec>(0, size-1, size);
 arma::uvec out = Rcpp::RcppArmadillo::sample(sequence, length, true, p);
 return out;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
// FORWARD FILTERING
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward Filtering - SingleStep
//[[Rcpp::export(name = "forward_filter")]]
Rcpp::List forward_filter(const arma::mat& Y, const arma::mat& G, const arma::mat& P, const arma::mat& V, const arma::mat& W, const arma::mat& m, const arma::mat& C, const double& nu, const arma::mat& Psi) {

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
Rcpp::List parallel_forward_filter(const arma::mat& Y, const arma::mat& G, const arma::mat& P, const arma::mat& D,
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
    } return results;
  }



///////////////////////////////////////////////////////////////////////////////////////////////////////////
// BACKWARD SAMPLING
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Backward Sampling - SingleStep
// [[Rcpp::export]]
arma::cube backward_sample(const arma::mat& G, const arma::mat& P, const arma::mat& V, const arma::mat& W, const List& ForwFilt, const arma::cube& ThetaSmp, const arma::cube& SigmaSmp,const int& t, const int& L){

  // Extract FF at time t
  List FF_t = ForwFilt(t);
  arma::mat C_inv = as<arma::mat>(FF_t["iC"]);
  arma::mat m = as<arma::mat>(FF_t["m"]);

  // Precomputing
  arma::mat G_trans = arma::trans(G);
  arma::mat W_inv = arma::inv(W);

  // Smoother posterior parameters
  arma::mat H_inv = C_inv + (G_trans * W_inv * G);
  arma::mat H = arma::inv(H_inv);

  // Set the environment
  Rcpp::Environment mniw = Rcpp::Environment::namespace_env("mniw");
  Rcpp::Function rMNorm_R = mniw["rMNorm"];
  arma::cube Theta_ss(ThetaSmp.n_rows,ThetaSmp.n_cols,L);
  Theta_ss.zeros();
  for (int l = 0; l < L; l++) {

    arma::mat Theta_l = ThetaSmp.slice(l);
    arma::mat Sigma_l = SigmaSmp.slice(l);
    arma::mat h = H * (C_inv * m + (G_trans * W_inv * Theta_l));
    Theta_ss.slice(l) = as<arma::mat>(rMNorm_R(Named("n", 1), Named("Lambda", h), Named("SigmaR", H), Named("SigmaC", Sigma_l)));

  }

  return Theta_ss;

}


// Backward Sampling - MultipleStep
// [[Rcpp::export]]
Rcpp::List backward_sample_T(const arma::mat& G, const arma::mat& P, const arma::mat& V, const arma::mat& W, const List& ForwFilt, const int& L){

  // Number of times
  int Tmax = ForwFilt.size()-1;

  // Set the environment
  Rcpp::Environment mniw = Rcpp::Environment::namespace_env("mniw");
  Rcpp::Function rMNIW_R = mniw["rMNIW"];

  // build containers
  List out_BS(Tmax+1);

  // LastStep - Filtered Posterior Information
  List out_FFtmax = ForwFilt(Tmax);
  arma::mat m_tmax = as<arma::mat>(out_FFtmax["m"]);
  arma::mat C_tmax = as<arma::mat>(out_FFtmax["C"]);
  double nu_tmax = as<double>(out_FFtmax["nu"]);
  arma::mat Psi_tmax = as<arma::mat>(out_FFtmax["Psi"]);
  List ThetaSigmaT = as<List>(rMNIW_R(Named("n", L), Named("Lambda", m_tmax), Named("Sigma", C_tmax), Named("nu", nu_tmax), Named("Psi", Psi_tmax)));
  arma::cube ThetaSmpT = as<arma::cube>(ThetaSigmaT["X"]);
  arma::cube SigmaSmpT = as<arma::cube>(ThetaSigmaT["V"]);
  arma::cube out_BStmax = backward_sample(G, P, V, W, ForwFilt, ThetaSmpT, SigmaSmpT, Tmax, L);
  out_BS(Tmax) = out_BStmax;

  // OtherSteps - Reverse For Loop
  for (int t = (Tmax-1); t >= 0; t--) {

    // Extract FF at time (t-1)
    List out_FF = ForwFilt(t);

    // Smoothed Sample
    arma::cube ThetaSmpBS = out_BS(t+1);
    out_BS(t) = backward_sample(G, P, V, W, ForwFilt, ThetaSmpBS, SigmaSmpT, t, L);

  }

  // Return List of Results
  return out_BS;

}


// Weighted Backward Sampling - SingleStep
// [[Rcpp::export]]
arma::cube weighted_backward_sample(const arma::mat& G, const arma::mat& D,
                                    const List& FF_t, const arma::cube& ThetaSmp, const arma::cube& SigmaSmp,
                                    const arma::mat& par_grid, const arma::vec& weights){

  // Gather information from inputs
  int n = D.n_rows;
  int p = G.n_rows - n;
  int J = par_grid.n_rows;
  int L = ThetaSmp.n_slices;
  // Rcpp::Rcout << "n: " << n << std::endl;
  // Rcpp::Rcout << "p: " << p << std::endl;
  // Rcpp::Rcout << "J: " << J << std::endl;
  // Rcpp::Rcout << "L: " << L << std::endl;

  // Sample the models indexes
  arma::uvec model_idx = sample_index(J, L, weights);
  arma::uvec uniqmod = arma::unique(model_idx);
  arma::uvec uniqind = find_unique(model_idx);
  // Rcpp::Rcout << "uniqmod: " << uniqmod << std::endl;
  // Rcpp::Rcout << "uniqind: " << uniqind << std::endl;
  int JJ = uniqmod.size();

  // Precomputing
  arma::mat G_trans = arma::trans(G);

  // Set the environment
  Rcpp::Environment mniw = Rcpp::Environment::namespace_env("mniw");
  Rcpp::Function rMNorm_R = mniw["rMNorm"];
  arma::cube Theta_ss(ThetaSmp.n_rows,ThetaSmp.n_cols,L);

  // Models loop
  List Models(JJ);
  for (int j = 0; j < JJ; j++) {

    // Sample the model
    arma::uword j_mod = uniqmod(j);
    List FF_j = FF_t(j_mod);
    arma::mat C_inv = as<arma::mat>(FF_j["iC"]);

    // Compute model specific objects
    double phi_j = par_grid(j, 1);
    arma::mat K_j = exp(-phi_j * D);
    arma::mat W_j = join_cols(
      join_rows(eye(p, p), zeros(p, n)),
      join_rows(zeros(n, p), K_j)
    );
    arma::mat W_inv = arma::inv(W_j);
    double tau_j = par_grid(j, 0);
    arma::mat V_j = ((1 - tau_j) / tau_j) * eye(n, n);

    // Smoother posterior parameters
    arma::mat H_inv = C_inv + (G_trans * W_inv * G);
    arma::mat H = arma::inv(H_inv);

    Models(j) = List::create(
      Named("m") = as<arma::mat>(FF_j["m"]),
      Named("C_inv") = C_inv,
      Named("H") = H,
      Named("W_inv") = W_inv);
  }

  // Rcpp::Rcout << "ok" << std::endl;

  // Samples loop
  Theta_ss.zeros();
  for (int l = 0; l < L; l++) {

    // Sample the model
    arma::uword l_mod = model_idx(l);
    arma::uvec lmod_idx = arma::find(uniqmod == l_mod);
    arma::uword lmod = lmod_idx(0);
    // Rcpp::Rcout << "l mod: " << l_mod << std::endl;
    // Rcpp::Rcout << "Models.size: " << Models.size() << std::endl;
    // Rcpp::Rcout << "l mod idx: " << lmod_idx << std::endl;
    // Rcpp::Rcout << "l mod: " << lmod << std::endl;
    List ModJ = Models(lmod);
    arma::mat m = ModJ["m"];
    arma::mat C_inv = ModJ["C_inv"];
    arma::mat W_inv = ModJ["W_inv"];
    arma::mat H = ModJ["H"];

    // Smoother posterior samples
    arma::mat Theta_l = ThetaSmp.slice(l);
    arma::mat Sigma_l = SigmaSmp.slice(l);
    arma::mat h = H * (C_inv * m + (G_trans * W_inv * Theta_l));
    Theta_ss.slice(l) = as<arma::mat>(rMNorm_R(Named("n", 1), Named("Lambda", h), Named("SigmaR", H), Named("SigmaC", Sigma_l)));

  }

  return Theta_ss;

}


// Weighted Backward Sampling - MultipleStep
// [[Rcpp::export]]
Rcpp::List weighted_backward_sample_T(const arma::mat& G, const arma::mat& D,
                                const List& ForwFilt, const int& L,
                                const arma::mat& par_grid, const arma::vec& weights) {

  // Gather Information from Inputs
  int J = par_grid.n_rows;
  int Tmax = ForwFilt.size()-1;

  // Set the Environment
  Rcpp::Environment mniw = Rcpp::Environment::namespace_env("mniw");
  Rcpp::Function rMNIW_R = mniw["rMNIW"];

  // Results Containers
  List out_BS(Tmax+1);

  // LastStep - Filtered Posterior Information
  List out_FFtmax = ForwFilt(Tmax);
  List FF_T = out_FFtmax["filtered_results"];

  // Sample the Models Indexes
  arma::uvec last_idx = sample_index(J, L, weights);
  arma::uvec uniqmod = arma::unique(last_idx);
  arma::uvec uniqind = find_unique(last_idx);
  int JJ = uniqmod.size();
  // Rcpp::Rcout << "last idx: " << last_idx << std::endl;
  // Rcpp::Rcout << "uniq mod: " << uniqmod << std::endl;
  // Rcpp::Rcout << "uniq ind: " << uniqind << std::endl;
  // Rcpp::Rcout << "JJ: " << JJ << std::endl;

  // Number of samples
  arma::uvec LLv = arma::hist(last_idx);
  arma::uvec LL_indices = arma::find(LLv > 0);
  arma::uvec LLi = LLv(LL_indices);

  arma::cube ThetaSmpT;
  arma::cube SigmaSmpT;
  for (int j = 0; j < JJ; j++) {

    // int j = 1;

    // Sample the Model
    arma::uword j_mod = uniqmod(j);
    List FF_j = FF_T(j_mod);

    // Rcpp::Rcout << "FF_j: " << FF_j << std::endl;

    arma::mat m_j = as<arma::mat>(FF_j["m"]);
    arma::mat C_j = as<arma::mat>(FF_j["C"]);
    double nu_j = as<double>(FF_j["nu"]);
    arma::mat Psi_j = as<arma::mat>(FF_j["Psi"]);

    // Number of Samples
    int LL = LLi(j);
    // Rcpp::Rcout << "LLv: " << LLv << std::endl;
    // Rcpp::Rcout << "LLv ind: " << LL_indices << std::endl;
    // Rcpp::Rcout << "LLi: " << LLi << std::endl;
    // Rcpp::Rcout << "uniq mod: " << uniqmod << std::endl;
    // Rcpp::Rcout << "j mod: " << j_mod << std::endl;
    // Rcpp::Rcout << "LL: " << LL << std::endl;

    // Posterior Ssample
    List ThetaSigmaT = as<List>(rMNIW_R(Named("n", LL), Named("Lambda", m_j), Named("Sigma", C_j), Named("nu", nu_j), Named("Psi", Psi_j)));

    arma::cube ThetaSmpT_j = as<arma::cube>(ThetaSigmaT["X"]);
    arma::cube SigmaSmpT_j = as<arma::cube>(ThetaSigmaT["V"]);

    // Save Results
    ThetaSmpT = join_slices(ThetaSmpT, ThetaSmpT_j);
    SigmaSmpT = join_slices(SigmaSmpT, SigmaSmpT_j);

  }

  // LastStep - Filtered Posterior Sample
  arma::uvec shuffle_indices = arma::randperm(L);
  // Rcpp::Rcout << "shuffle ind: " << shuffle_indices << std::endl;
  arma::cube ThetaSmpTmax = ThetaSmpT.slices(shuffle_indices);
  arma::cube SigmaSmpTmax = SigmaSmpT.slices(shuffle_indices);
  arma::cube out_BStmax = weighted_backward_sample(G, D, FF_T, ThetaSmpTmax, SigmaSmpTmax, par_grid, weights);
  out_BS(Tmax) = out_BStmax;

  // OtherSteps - Reverse For Loop
  for (int t = (Tmax-1); t >= 0; t--) {

    // Extract FF at time (t-1)
    List ForwFilt_t = ForwFilt(t);
    List out_FF = ForwFilt_t["filtered_results"];

    // Smoothed Sample
    arma::cube ThetaSmpBS = out_BS(t+1);
    out_BS(t) = weighted_backward_sample(G, D, out_FF, ThetaSmpBS, SigmaSmpT, par_grid, weights);

  }

  // Return List of Results
  return out_BS;

}


// Weighted Backward Sampling - MultipleStep
// [[Rcpp::export]]
Rcpp::List weighted_backward_sample_T3(const arma::cube& G, const arma::mat& D,
                                 const List& ForwFilt, const int& L,
                                 const arma::mat& par_grid, const arma::vec& weights) {

  // Gather Information from Inputs
  int J = par_grid.n_rows;
  int Tmax = ForwFilt.size()-1;

  // Set the Environment
  Rcpp::Environment mniw = Rcpp::Environment::namespace_env("mniw");
  Rcpp::Function rMNIW_R = mniw["rMNIW"];

  // Results Containers
  List out_BS(Tmax+1);

  // LastStep - Filtered Posterior Information
  List out_FFtmax = ForwFilt(Tmax);
  List FF_T = out_FFtmax["filtered_results"];

  // Sample the Models Indexes
  arma::uvec last_idx = sample_index(J, L, weights);
  arma::uvec uniqmod = arma::unique(last_idx);
  arma::uvec uniqind = find_unique(last_idx);
  int JJ = uniqmod.size();
  // Rcpp::Rcout << "last idx: " << last_idx << std::endl;
  // Rcpp::Rcout << "uniq mod: " << uniqmod << std::endl;
  // Rcpp::Rcout << "uniq ind: " << uniqind << std::endl;
  // Rcpp::Rcout << "JJ: " << JJ << std::endl;

  // Number of samples
  arma::uvec LLv = arma::hist(last_idx);
  arma::uvec LL_indices = arma::find(LLv > 0);
  arma::uvec LLi = LLv(LL_indices);

  arma::cube ThetaSmpT;
  arma::cube SigmaSmpT;
  for (int j = 0; j < JJ; j++) {

    // int j = 1;

    // Sample the Model
    arma::uword j_mod = uniqmod(j);
    List FF_j = FF_T(j_mod);

    // Rcpp::Rcout << "FF_j: " << FF_j << std::endl;

    arma::mat m_j = as<arma::mat>(FF_j["m"]);
    arma::mat C_j = as<arma::mat>(FF_j["C"]);
    double nu_j = as<double>(FF_j["nu"]);
    arma::mat Psi_j = as<arma::mat>(FF_j["Psi"]);

    // Number of Samples
    int LL = LLi(j);
    // Rcpp::Rcout << "LLv: " << LLv << std::endl;
    // Rcpp::Rcout << "LLv ind: " << LL_indices << std::endl;
    // Rcpp::Rcout << "LLi: " << LLi << std::endl;
    // Rcpp::Rcout << "uniq mod: " << uniqmod << std::endl;
    // Rcpp::Rcout << "j mod: " << j_mod << std::endl;
    // Rcpp::Rcout << "LL: " << LL << std::endl;

    // Posterior Ssample
    List ThetaSigmaT = as<List>(rMNIW_R(Named("n", LL), Named("Lambda", m_j), Named("Sigma", C_j), Named("nu", nu_j), Named("Psi", Psi_j)));

    arma::cube ThetaSmpT_j = as<arma::cube>(ThetaSigmaT["X"]);
    arma::cube SigmaSmpT_j = as<arma::cube>(ThetaSigmaT["V"]);

    // Save Results
    ThetaSmpT = join_slices(ThetaSmpT, ThetaSmpT_j);
    SigmaSmpT = join_slices(SigmaSmpT, SigmaSmpT_j);

  }

  // LastStep - Filtered Posterior Sample
  arma::uvec shuffle_indices = arma::randperm(L);
  // Rcpp::Rcout << "shuffle ind: " << shuffle_indices << std::endl;
  arma::cube ThetaSmpTmax = ThetaSmpT.slices(shuffle_indices);
  arma::cube SigmaSmpTmax = SigmaSmpT.slices(shuffle_indices);
  arma::mat GT = G.slice(Tmax);
  arma::cube out_BStmax = weighted_backward_sample(GT, D, FF_T, ThetaSmpTmax, SigmaSmpTmax, par_grid, weights);
  out_BS(Tmax) = out_BStmax;

  // OtherSteps - Reverse For Loop
  for (int t = (Tmax-1); t >= 0; t--) {

    // Extract FF at time (t-1)
    List ForwFilt_t = ForwFilt(t);
    List out_FF = ForwFilt_t["filtered_results"];

    // Smoothed Sample
    arma::cube ThetaSmpBS = out_BS(t+1);
    arma::mat Gt = G.slice(t);
    out_BS(t) = weighted_backward_sample(Gt, D, out_FF, ThetaSmpBS, SigmaSmpT, par_grid, weights);

  }

  // Return List of Results
  return out_BS;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
// MATRIX T EVALUATIONS
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Multivariate Gamma function
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

// Projects a vector onto the probability simplex
std::vector<double> project_simplex(std::vector<double> v) {
  int n = v.size();
  std::vector<double> u(v);
  std::sort(u.begin(), u.end(), std::greater<double>());
  double sum = 0.0;
  double rho = 0.0;

  for (int j = 0; j < n; ++j) {
    sum += u[j];
    double t = (sum - 1.0) / (j + 1);
    if (u[j] > t) {
      rho = j + 1;
    } else {
      break;
    }
  }

  double tau = (std::accumulate(u.begin(), u.begin() + rho, 0.0) - 1.0) / rho;
  std::vector<double> w(n);
  for (int i = 0; i < n; ++i)
    w[i] = std::max(v[i] - tau, 0.0);
  return w;
}

// Gradient Descent - Projected onto simplex
// [[Rcpp::export]]
NumericVector optimize_weights_proj(NumericMatrix scores, double lr = 0.05, int max_iter = 500) {
  int n = scores.nrow(), J = scores.ncol();
  std::vector<double> w(J, 1.0 / J);  // uniform init

  for (int iter = 0; iter < max_iter; ++iter) {
    // Gradient of the mean log-sum
    std::vector<double> grad(J, 0.0);

    for (int i = 0; i < n; ++i) {
      double dot = 0.0;
      for (int j = 0; j < J; ++j)
        dot += scores(i, j) * w[j];

      for (int j = 0; j < J; ++j)
        grad[j] += scores(i, j) / dot;
    }

    for (int j = 0; j < J; ++j)
      grad[j] /= n;

    // Gradient ascent step
    for (int j = 0; j < J; ++j)
      w[j] += lr * grad[j];

    // Project back onto the simplex
    w = project_simplex(w);
  }

  return NumericVector(w.begin(), w.end());
}


// Adam Optimization - Porjected onto simplex
// [[Rcpp::export]]
NumericVector optimize_weights_adam(NumericMatrix scores, double lr = 0.05, int max_iter = 1000) {
  int n = scores.nrow(), J = scores.ncol();
  std::vector<double> w(J, 1.0 / J);  // uniform init
  std::vector<double> m(J, 0.0), v(J, 0.0);
  double beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8;

  for (int t = 1; t <= max_iter; ++t) {
    std::vector<double> grad(J, 0.0);

    for (int i = 0; i < n; ++i) {
      double dot = 0.0;
      for (int j = 0; j < J; ++j)
        dot += scores(i, j) * w[j];
      for (int j = 0; j < J; ++j)
        grad[j] += scores(i, j) / dot;
    }
    for (int j = 0; j < J; ++j)
      grad[j] /= n;

    // Adam updates
    for (int j = 0; j < J; ++j) {
      m[j] = beta1 * m[j] + (1 - beta1) * grad[j];
      v[j] = beta2 * v[j] + (1 - beta2) * grad[j] * grad[j];

      double m_hat = m[j] / (1 - std::pow(beta1, t));
      double v_hat = v[j] / (1 - std::pow(beta2, t));

      w[j] += lr * m_hat / (std::sqrt(v_hat) + epsilon);
    }

    // Project onto simplex
    w = project_simplex(w);
  }

  return NumericVector(w.begin(), w.end());
}


// [[Rcpp::export]]
NumericMatrix compute_Wt_cpp(Rcpp::List density_list,
                             int n, // number of observations
                             int t, // number of time steps (T)
                             double lr = 0.05,
                             int max_iter = 500,
                             int n_threads = 0) {
  Rcout << "Computing Dynamic Bayesian Predictive Stacking Weights (C++) ...\n";

  if (n_threads > 0)
    omp_set_num_threads(n_threads);

  int J = as<NumericMatrix>(density_list[0]).ncol();  // number of models
  NumericMatrix Wi(n, J);  // result matrix

#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    NumericMatrix epd_i(t - 1, J);
    for (int a = 1; a < t; ++a) {  // from time 2 to T
      NumericMatrix mat = density_list[a - 1];  // density_list[0] is time 2
      for (int j = 0; j < J; ++j)
        epd_i(a - 1, j) = std::exp(mat(i, j));  // exp(t(epd_i))
    }

    Rcpp::NumericVector w = optimize_weights_proj(epd_i, lr, max_iter);
    for (int j = 0; j < J; ++j)
      Wi(i, j) = w[j];
  }

  return Wi;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
// UNIFIED FUNCTION
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Unified Parallel Function
// [[Rcpp::export]]
Rcpp::List unified_parallel_function(const arma::mat& Y, const arma::mat& G, const arma::mat& P, const arma::mat& D,
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


// Parallel BPS-Forward Filtering - Multiple Step
// [[Rcpp::export]]
Rcpp::List spFF(const arma::cube& Y, const arma::mat& G, const arma::mat& P, const arma::mat& D,
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
    // arma::mat Yt = Y.slice(t);
    // List out_old = out_FF(t-1);
    // arma::mat m_old = as<arma::mat>(out_old["m"]);
    // arma::mat C_old = as<arma::mat>(out_old["C"]);
    // double nu_old = as<double>(out_old["nu"]);
    // arma::mat Psi_old = as<arma::mat>(out_old["Psi"]);
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


// Parallel BPS-Forward Filtering - Multiple Step
// [[Rcpp::export]]
Rcpp::List spFF3(const arma::cube& Y, const arma::cube& G, const arma::cube& P, const arma::mat& D,
           const arma::mat& par_grid, List const& prior,
           int num_threads = 1) {

  // Rcpp::Rcout << "Checkpoint 1" << std::endl;

  // Gather information from inputs
  int n = D.n_rows;
  int p = G.n_rows - n;
  int T = Y.n_slices;
  int J = par_grid.n_rows;

  // build containers
  List out_FF(T);

  // Rcpp::Rcout << "Checkpoint 2" << std::endl;

  // FirstStep - Prior Information
  // arma::mat Y1 = Y.slice(0);
  // arma::mat m0 = as<arma::mat>(prior["m"]);
  // arma::mat C0 = as<arma::mat>(prior["C"]);
  // double nu0 = as<double>(prior["nu"]);
  // arma::mat Psi0 = as<arma::mat>(prior["Psi"]);
  // List out_prior = forward_filter(Y1, G, P, V, W, m0, C0, nu0, Psi0);
  arma::mat Y1 = Y.slice(0);
  arma::mat G1 = G.slice(0);
  arma::mat P1 = P.slice(0);
  std::vector<List> priors = std::vector<List>(J, prior);
  List prior_list = Rcpp::wrap(priors);

  // Rcpp::Rcout << "Checkpoint 3" << std::endl;
  List out_prior = parallel_forward_filter(Y1, G1, P1, D, par_grid, prior_list, J, p, n, num_threads);
  out_FF(0) = List::create(Named("filtered_results") = out_prior);

  // Rcpp::Rcout << "Checkpoint 4" << std::endl;
  // OtherSteps - For Loop
  for (int t = 1; t < T; t++) {

    // Extract Old Information
    // arma::mat Yt = Y.slice(t);
    // List out_old = out_FF(t-1);
    // arma::mat m_old = as<arma::mat>(out_old["m"]);
    // arma::mat C_old = as<arma::mat>(out_old["C"]);
    // double nu_old = as<double>(out_old["nu"]);
    // arma::mat Psi_old = as<arma::mat>(out_old["Psi"]);
    arma::mat Yt = Y.slice(t);
    arma::mat Gt = G.slice(t);
    arma::mat Pt = P.slice(t);
    List out_prec = out_FF(t-1);
    List out_old = out_prec["filtered_results"];

    // Compute New Information
    List out_new = unified_parallel_function(Yt, Gt, Pt, D, par_grid, out_old, J, p, n, num_threads);
    out_FF(t) = out_new;

  }

  // Rcpp::Rcout << "Checkpoint 5" << std::endl;
  // Return List of Results
  return out_FF;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
// TEMPORAL FORECASTING
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// 1-Step Ahead Posterior Predictive Samples
// [[Rcpp::export]]
Rcpp::List predict_ffbs(const arma::mat& G, const arma::mat& P, const arma::mat& V, const arma::mat& W,
                 const arma::mat& a, const arma::mat& R, const double& nu, const arma::mat& Psi,
                 const int& L) {

 // Pre-computing
 arma::mat Gt = arma::trans(G);
 arma::mat Pt = arma::trans(P);
 // Rcpp::Rcout << "ok! " << std::endl;

 // Predictive parameters
 arma::mat a_new = G * a;
 arma::mat R_new = G * R * Gt + W;
 // Rcpp::Rcout << "ok! " << std::endl;

 // 1-Step ahead predictive
 arma::mat f_pred = P * a_new;
 arma::mat Q_pred = P * R_new * Pt + V;
 // Rcpp::Rcout << "ok! " << std::endl;

 // Set the environment
 Rcpp::Environment mniw = Rcpp::Environment::namespace_env("mniw");
 Rcpp::Function rMT_R = mniw["rMT"];

 // Predictive samples
 arma::cube Y_pred = as<arma::cube>(rMT_R(Named("n", L), Named("Lambda", f_pred), Named("SigmaR", Q_pred), Named("SigmaC", Psi), Named("nu", nu)));

 // Return list
 return List::create(Named("Y") = Y_pred,
                     Named("a") = a_new,
                     Named("R") = R_new);

}


// 1-Step Ahead Forecast Samples
// [[Rcpp::export]]
Rcpp::List temporal_forecast(const arma::mat& G, const arma::mat& P,
                      const arma::mat& D, const arma::mat& par_grid,
                      const List& ForwFilt, const arma::mat& weights,
                      const int& horiz, const int& L) {

 // Gather data dimensions
 int n = D.n_rows;
 int p = G.n_rows - n;
 int J = par_grid.n_rows;
 int tmax = weights.n_cols + 1;
 int th = tmax + horiz;

 // Result container
 List pred(th);
 List predictions(th);

 // Rcpp::Rcout << "ok 1 " << std::endl;

 for (int t = 0; t < th; t++) {

   // Select the weights
   arma::uword w_idx = std::min(t, tmax-1);
   // Rcpp::Rcout << "w_idx: " << w_idx << std::endl;
   // Rcpp::Rcout << "w ncol: " << weights.n_cols << std::endl;
   arma::vec weights_t;
   if(t < 1) {
     arma::vec weights_t = arma::vec(J, arma::fill::ones) / J;
   } else {
     arma::vec weights_t = weights.col(w_idx-1);
     // Rcpp::Rcout << "weights: " << weights.col(w_idx-1) << std::endl;
   }

   // Sample the model
   arma::uvec model_idx = sample_index(J, 1, weights_t);
   arma::uword j = model_idx(0);

   // Compute model specific objects
   double phi_j = par_grid(j, 1);
   arma::mat K_j = exp(-phi_j * D);
   arma::mat W_j = join_cols(
     join_rows(eye(p, p), zeros(p, n)),
     join_rows(zeros(n, p), K_j)
   );
   double tau_j = par_grid(j, 0);
   arma::mat V_j = ((1 - tau_j) / tau_j) * eye(n, n);

   // Rcpp::Rcout << "ok 2 " << std::endl;

   // Predictive samples
   if (t < tmax) {
     // Extract FF results
     List FF_t = ForwFilt(t);
     List FF_J = FF_t["filtered_results"];
     List FF_j = FF_J(j);
     List pred_t = predict_ffbs(G, P, V_j, W_j, FF_j["m"], FF_j["C"], FF_j["nu"], FF_j["Psi"], L);
     pred[t] = pred_t;
     predictions[t] = pred_t["Y"];
     // Rcpp::Rcout << "ok 3 " << std::endl;
   } else {
     // Extract Previous results
     List pred_prev = pred[t-1];
     List FF_t = ForwFilt(tmax-1);
     List FF_J = FF_t["filtered_results"];
     List FF_j = FF_J(j);
     List pred_t = predict_ffbs(G, P, V_j, W_j, pred_prev["a"], pred_prev["R"], FF_j["nu"], FF_j["Psi"], L);
     pred[t] = pred_t;
     predictions[t] = pred_t["Y"];
     // Rcpp::Rcout << "ok 4 " << std::endl;
   }

 }

 // List return
 return predictions;

}


// 1-Step Ahead Forecast Samples
// [[Rcpp::export]]
Rcpp::List temporal_forecast3(const arma::cube& G, const arma::cube& P,
                       const arma::mat& D, const arma::mat& par_grid,
                       const List& ForwFilt, const arma::mat& weights,
                       const int& horiz, const int& L) {

 // Gather data dimensions
 int n = D.n_rows;
 int p = G.n_rows - n;
 int J = par_grid.n_rows;
 int tmax = weights.n_cols + 1;
 int th = tmax + horiz;

 // Result container
 List pred(th);
 List predictions(th);

 // Rcpp::Rcout << "ok 1 " << std::endl;

 for (int t = 0; t < th; t++) {

   // Select the weights
   arma::uword w_idx = std::min(t, tmax-1);
   // Rcpp::Rcout << "w_idx: " << w_idx << std::endl;
   // Rcpp::Rcout << "w ncol: " << weights.n_cols << std::endl;
   arma::vec weights_t;
   if(t < 1) {
     arma::vec weights_t = arma::vec(J, arma::fill::ones) / J;
   } else {
     arma::vec weights_t = weights.col(w_idx-1);
     // Rcpp::Rcout << "weights: " << weights.col(w_idx-1) << std::endl;
   }

   // Sample the model
   arma::uvec model_idx = sample_index(J, 1, weights_t);
   arma::uword j = model_idx(0);

   // Compute model specific objects
   double phi_j = par_grid(j, 1);
   arma::mat K_j = exp(-phi_j * D);
   arma::mat W_j = join_cols(
     join_rows(eye(p, p), zeros(p, n)),
     join_rows(zeros(n, p), K_j)
   );
   double tau_j = par_grid(j, 0);
   arma::mat V_j = ((1 - tau_j) / tau_j) * eye(n, n);

   // Rcpp::Rcout << "ok 2 " << std::endl;

   // Predictive samples
   if (t < tmax) {
     // Extract FF results
     List FF_t = ForwFilt(t);
     List FF_J = FF_t["filtered_results"];
     List FF_j = FF_J(j);
     arma::mat Gt = G.slice(t);
     arma::mat Pt = P.slice(t);
     List pred_t = predict_ffbs(Gt, Pt, V_j, W_j, FF_j["m"], FF_j["C"], FF_j["nu"], FF_j["Psi"], L);
     pred[t] = pred_t;
     predictions[t] = pred_t["Y"];
     // Rcpp::Rcout << "ok 3 " << std::endl;
   } else {
     // Extract Previous results
     List pred_prev = pred[t-1];
     List FF_t = ForwFilt(tmax-1);
     List FF_J = FF_t["filtered_results"];
     List FF_j = FF_J(j);
     arma::mat Gt = G.slice(t-1);
     arma::mat Pt = P.slice(t-1);
     List pred_t = predict_ffbs(Gt, Pt, V_j, W_j, pred_prev["a"], pred_prev["R"], FF_j["nu"], FF_j["Psi"], L);
     pred[t] = pred_t;
     predictions[t] = pred_t["Y"];
     // Rcpp::Rcout << "ok 4 " << std::endl;
   }

 }

 // List return
 return predictions;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
// SPATIAL INTERPOLATION
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Spatial Interpolation
// [[Rcpp::export]]
Rcpp::List spatial_interpolation(const arma::mat& G, const arma::mat& P, const arma::mat& Xu,
                          const arma::mat& D_s, const arma::mat& D_u, const arma::mat& D_us,
                          const arma::mat& par_grid,
                          const List& ForwFilt, const arma::mat& weights,
                          const int& t, const int& L) {

 // Gather data dimensions
 int n = D_s.n_rows;
 int u = D_u.n_rows;
 int p = G.n_rows - n;
 int J = par_grid.n_rows;

 // Result container
 List predictions(L);
 //  TBD (marginal predictive sampling)
 // List Y_pred(L);
 // List O_pred(L);

 // Select the weights
 arma::vec w0 = arma::vec(J, arma::fill::ones) / J;
 arma::mat weightst = arma::join_horiz(w0, weights);
 arma::vec weights_t = weightst.col(t-1);
 // Rcpp::Rcout << weightst.n_cols << std::endl;
 // Rcpp::Rcout << t-1 << std::endl;

 // arma::vec weights_t;
 // if((t-1) < 1) {
 //   arma::vec weights_t = arma::vec(J, arma::fill::ones) / J;
 // } else {
 //   arma::vec weights_t = weights.col(t-1);
 // }

 // Sample the models indexes
 arma::uvec model_idx = sample_index(J, L, weights_t);
 arma::uvec uniqmod = arma::unique(model_idx);
 arma::uvec uniqind = find_unique(model_idx);
 int JJ = uniqmod.size();

 // Set the environment
 Rcpp::Environment mniw = Rcpp::Environment::namespace_env("mniw");
 Rcpp::Function rMT_R = mniw["rMT"];

 // Rcpp::Rcout << "ok 1 " << std::endl;
 // Models loop
 List Models(JJ);
 for (int j = 0; j < JJ; j++) {
   // Sample the model
   arma::uword j_mod = uniqmod(j);

   // Compute model specific objects
   double tau_j = par_grid(j_mod, 0);
   arma::mat Vtilde = ((1 - tau_j) / tau_j) * eye(u, u);

   double phi_j = par_grid(j_mod, 1);
   arma::mat Kss = exp(-phi_j * D_s);
   arma::mat Kuu = exp(-phi_j * D_u);
   arma::mat Kus = exp(-phi_j * D_us);
   arma::mat iKss = arma::inv(Kss);

   arma::mat Mtilde = Kus * iKss;
   arma::mat Wtilde = Kuu - (Mtilde * arma::trans(Kus));

   List FF_t = ForwFilt(t-1);
   List FF_J = FF_t["filtered_results"];
   List FF_j = FF_J(j_mod);
   arma::mat C_t = FF_j["C"];
   arma::mat m_t = FF_j["m"];
   arma::mat Psi_t = FF_j["Psi"];
   double nu_t = FF_j["nu"];

   arma::mat zero_up = arma::zeros(u, p);
   arma::mat chi_t = arma::join_horiz(arma::join_vert(Xu, zero_up), arma::join_vert(Mtilde, Mtilde));
   arma::mat N_t = arma::join_horiz(arma::join_vert(Vtilde + Wtilde, Wtilde), arma::join_vert(Wtilde, Wtilde));

   arma::mat mu_t = chi_t * m_t;
   arma::mat E_t = chi_t * C_t * arma::trans(chi_t) + N_t;

   Models(j) = List::create(
     Named("mu") = mu_t,
     Named("E") = E_t,
     Named("nu") = nu_t,
     Named("Psi") = Psi_t);
 }

 // Rcpp::Rcout << "ok 1 " << std::endl;
 // Samples loop
 for (int l = 0; l < L; l++) {

   // Sample the model
   arma::uword l_mod = model_idx(l);
   arma::uvec lmod_idx = arma::find(uniqmod == l_mod);
   arma::uword lmod = lmod_idx(0);
   List ModJ = Models(lmod);

   arma::mat mu = as<arma::mat>(ModJ["mu"]);
   arma::mat E = as<arma::mat>(ModJ["E"]);
   arma::mat Psi = as<arma::mat>(ModJ["Psi"]);
   double nu = as<double>(ModJ["nu"]);

   // Predictive sample
   predictions[l] = as<arma::mat>(rMT_R(Named("n", 1), Named("Lambda", mu), Named("SigmaR", E), Named("SigmaC", Psi), Named("nu", nu)));

 }

 // List return
 return predictions;

}
