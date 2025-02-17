#include <RcppArmadillo.h>
#include <osqp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// // Function to solve the convex optimization problem
// // [[Rcpp::export]]
// arma::vec conv_opt_rcpp(const arma::mat& scores) {
//   int n = scores.n_rows;
//   int K = scores.n_cols;
//
//   // OSQP data
//   arma::sp_mat P(K, K);  // Quadratic term (zero in this case)
//   arma::vec q(K, fill::zeros);  // Linear term
//   arma::sp_mat A(1 + K, K);  // Constraint matrix
//   arma::vec l(1 + K);  // Lower bounds
//   arma::vec u(1 + K);  // Upper bounds
//
//   // Initialize P and q
//   P.zeros();
//   q.zeros();
//
//   // Constraints: weights >= 0 and sum(weights) == 1
//   A.zeros();
//   A(0, 0) = 1.0;  // Sum constraint
//   for (int i = 0; i < K; ++i) {
//     A(1 + i, i) = 1.0;  // Non-negativity constraints
//   }
//
//   l(0) = 1.0;  // Sum of weights == 1
//   u(0) = 1.0;
//   for (int i = 1; i <= K; ++i) {
//     l(i) = 0.0;  // Weights >= 0
//     u(i) = OSQP_INFTY;
//   }
//
//   // OSQP settings
//   OSQPSettings settings;
//   osqp_set_default_settings(&settings);
//   settings.verbose = false;
//
//   // OSQP workspace
//   OSQPWorkspace* work = osqp_setup(P.n_rows, A.n_rows, P.values, P.row_indices, P.col_ptrs,
//                                    q.memptr(), A.values, A.row_indices, A.col_ptrs,
//                                    l.memptr(), u.memptr(), &settings);
//
//   // Solve the problem
//   osqp_solve(work);
//
//   // Extract the solution
//   arma::vec weights(K);
//   for (int i = 0; i < K; ++i) {
//     weights(i) = work->solution->x[i];
//   }
//
//   // Clean up
//   osqp_cleanup(work);
//
//   return weights;
// }
