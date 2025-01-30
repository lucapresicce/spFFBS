// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// forward_filter
List forward_filter(const arma::mat& Y, const arma::mat& G, const arma::mat& P, const arma::mat& V, const arma::mat& W, const arma::mat& m, const arma::mat& C, const double& nu, const arma::mat& Psi);
RcppExport SEXP _spFFBS_forward_filter(SEXP YSEXP, SEXP GSEXP, SEXP PSEXP, SEXP VSEXP, SEXP WSEXP, SEXP mSEXP, SEXP CSEXP, SEXP nuSEXP, SEXP PsiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const double& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Psi(PsiSEXP);
    rcpp_result_gen = Rcpp::wrap(forward_filter(Y, G, P, V, W, m, C, nu, Psi));
    return rcpp_result_gen;
END_RCPP
}
// forward_filter_T
List forward_filter_T(const arma::cube& Y, const arma::mat& G, const arma::mat& P, const arma::mat& V, const arma::mat& W, List const& prior);
RcppExport SEXP _spFFBS_forward_filter_T(SEXP YSEXP, SEXP GSEXP, SEXP PSEXP, SEXP VSEXP, SEXP WSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< List const& >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(forward_filter_T(Y, G, P, V, W, prior));
    return rcpp_result_gen;
END_RCPP
}
// backward_sample
arma::cube backward_sample(const arma::mat& G, const arma::mat& P, const arma::mat& V, const arma::mat& W, const List& ForwFilt, const arma::cube& ThetaSmp, const arma::cube& SigmaSmp, const int& t, const int& L);
RcppExport SEXP _spFFBS_backward_sample(SEXP GSEXP, SEXP PSEXP, SEXP VSEXP, SEXP WSEXP, SEXP ForwFiltSEXP, SEXP ThetaSmpSEXP, SEXP SigmaSmpSEXP, SEXP tSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const List& >::type ForwFilt(ForwFiltSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type ThetaSmp(ThetaSmpSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type SigmaSmp(SigmaSmpSEXP);
    Rcpp::traits::input_parameter< const int& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int& >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(backward_sample(G, P, V, W, ForwFilt, ThetaSmp, SigmaSmp, t, L));
    return rcpp_result_gen;
END_RCPP
}
// backward_sample_T
List backward_sample_T(const arma::mat& G, const arma::mat& P, const arma::mat& V, const arma::mat& W, const List& ForwFilt, const int& L);
RcppExport SEXP _spFFBS_backward_sample_T(SEXP GSEXP, SEXP PSEXP, SEXP VSEXP, SEXP WSEXP, SEXP ForwFiltSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const List& >::type ForwFilt(ForwFiltSEXP);
    Rcpp::traits::input_parameter< const int& >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(backward_sample_T(G, P, V, W, ForwFilt, L));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spFFBS_forward_filter", (DL_FUNC) &_spFFBS_forward_filter, 9},
    {"_spFFBS_forward_filter_T", (DL_FUNC) &_spFFBS_forward_filter_T, 6},
    {"_spFFBS_backward_sample", (DL_FUNC) &_spFFBS_backward_sample, 9},
    {"_spFFBS_backward_sample_T", (DL_FUNC) &_spFFBS_backward_sample_T, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_spFFBS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
