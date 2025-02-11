#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

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


// Forward Filtering - MultipleStep
//[[Rcpp::export]]
List forward_filter_T(const arma::cube& Y, const arma::mat& G, const arma::mat& P, const arma::mat& V, const arma::mat& W, List const& prior){

   // Define Number of Temporal Slices
   int T = Y.n_slices;

   // build containers
   List out_FF(T);

   // FirstStep - Prior Information
   arma::mat Y1 = Y.slice(0);
   arma::mat m0 = as<arma::mat>(prior["m"]);
   arma::mat C0 = as<arma::mat>(prior["C"]);
   double nu0 = as<double>(prior["nu"]);
   arma::mat Psi0 = as<arma::mat>(prior["Psi"]);
   List out_prior = forward_filter(Y1, G, P, V, W, m0, C0, nu0, Psi0);
   out_FF(0) = out_prior;

   // OtherSteps - For Loop
   for (int t = 1; t < T; t++) {

     // Extract Old Information
     arma::mat Yt = Y.slice(t);
     List out_old = out_FF(t-1);
     arma::mat m_old = as<arma::mat>(out_old["m"]);
     arma::mat C_old = as<arma::mat>(out_old["C"]);
     double nu_old = as<double>(out_old["nu"]);
     arma::mat Psi_old = as<arma::mat>(out_old["Psi"]);

     // Compute New Information
     List out_new = forward_filter(Yt, G, P, V, W, m_old, C_old, nu_old, Psi_old);
     out_FF(t) = out_new;

   }

   // Return List of Results
   return out_FF;

 }


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// BACKWARD SAMPLING
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Backward Sampling - SingleStep
//[[Rcpp::export]]
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
//[[Rcpp::export]]
List backward_sample_T(const arma::mat& G, const arma::mat& P, const arma::mat& V, const arma::mat& W, const List& ForwFilt, const int& L){

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
//[[Rcpp::export]]
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


//' Weighted Backward Sampling - MultipleStep
//'
//' @export
//[[Rcpp::export]]
List spBS(const arma::mat& G, const arma::mat& D,
          const List& ForwFilt, const int& L,
          const arma::mat& par_grid, const arma::vec& weights) {

   // Gather Inforarma::mation from Inputs
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
