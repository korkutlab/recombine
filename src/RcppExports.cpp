// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_outliers_from_dist
List get_outliers_from_dist(NumericMatrix dist_, int loop_k, double loop_lambda, double loop_threshold, bool outlier_on);
RcppExport SEXP _recombine_get_outliers_from_dist(SEXP dist_SEXP, SEXP loop_kSEXP, SEXP loop_lambdaSEXP, SEXP loop_thresholdSEXP, SEXP outlier_onSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dist_(dist_SEXP);
    Rcpp::traits::input_parameter< int >::type loop_k(loop_kSEXP);
    Rcpp::traits::input_parameter< double >::type loop_lambda(loop_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type loop_threshold(loop_thresholdSEXP);
    Rcpp::traits::input_parameter< bool >::type outlier_on(outlier_onSEXP);
    rcpp_result_gen = Rcpp::wrap(get_outliers_from_dist(dist_, loop_k, loop_lambda, loop_threshold, outlier_on));
    return rcpp_result_gen;
END_RCPP
}
// get_outliers_from_ds
List get_outliers_from_ds(NumericMatrix ds_, int loop_k, double loop_lambda, double loop_threshold, bool outlier_on);
RcppExport SEXP _recombine_get_outliers_from_ds(SEXP ds_SEXP, SEXP loop_kSEXP, SEXP loop_lambdaSEXP, SEXP loop_thresholdSEXP, SEXP outlier_onSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ds_(ds_SEXP);
    Rcpp::traits::input_parameter< int >::type loop_k(loop_kSEXP);
    Rcpp::traits::input_parameter< double >::type loop_lambda(loop_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type loop_threshold(loop_thresholdSEXP);
    Rcpp::traits::input_parameter< bool >::type outlier_on(outlier_onSEXP);
    rcpp_result_gen = Rcpp::wrap(get_outliers_from_ds(ds_, loop_k, loop_lambda, loop_threshold, outlier_on));
    return rcpp_result_gen;
END_RCPP
}
// SHC_lasso_getuw
List SHC_lasso_getuw(NumericMatrix ds_, NumericVector wbounds_, int max_iter, bool init_random, bool silent, bool warm_start);
RcppExport SEXP _recombine_SHC_lasso_getuw(SEXP ds_SEXP, SEXP wbounds_SEXP, SEXP max_iterSEXP, SEXP init_randomSEXP, SEXP silentSEXP, SEXP warm_startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ds_(ds_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wbounds_(wbounds_SEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type init_random(init_randomSEXP);
    Rcpp::traits::input_parameter< bool >::type silent(silentSEXP);
    Rcpp::traits::input_parameter< bool >::type warm_start(warm_startSEXP);
    rcpp_result_gen = Rcpp::wrap(SHC_lasso_getuw(ds_, wbounds_, max_iter, init_random, silent, warm_start));
    return rcpp_result_gen;
END_RCPP
}
// SHC_lasso_lagrange_getuw
List SHC_lasso_lagrange_getuw(NumericMatrix ds_, NumericVector lambdas_, int max_iter, bool init_random, bool silent, bool warm_start);
RcppExport SEXP _recombine_SHC_lasso_lagrange_getuw(SEXP ds_SEXP, SEXP lambdas_SEXP, SEXP max_iterSEXP, SEXP init_randomSEXP, SEXP silentSEXP, SEXP warm_startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ds_(ds_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambdas_(lambdas_SEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type init_random(init_randomSEXP);
    Rcpp::traits::input_parameter< bool >::type silent(silentSEXP);
    Rcpp::traits::input_parameter< bool >::type warm_start(warm_startSEXP);
    rcpp_result_gen = Rcpp::wrap(SHC_lasso_lagrange_getuw(ds_, lambdas_, max_iter, init_random, silent, warm_start));
    return rcpp_result_gen;
END_RCPP
}
// SHC_FL_getuw
List SHC_FL_getuw(NumericMatrix ds_, NumericVector lambda1s_, NumericVector lambda2s_, int max_iter, bool init_random, bool silent);
RcppExport SEXP _recombine_SHC_FL_getuw(SEXP ds_SEXP, SEXP lambda1s_SEXP, SEXP lambda2s_SEXP, SEXP max_iterSEXP, SEXP init_randomSEXP, SEXP silentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ds_(ds_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda1s_(lambda1s_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda2s_(lambda2s_SEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type init_random(init_randomSEXP);
    Rcpp::traits::input_parameter< bool >::type silent(silentSEXP);
    rcpp_result_gen = Rcpp::wrap(SHC_FL_getuw(ds_, lambda1s_, lambda2s_, max_iter, init_random, silent));
    return rcpp_result_gen;
END_RCPP
}
// SHC_SSL_getuw
List SHC_SSL_getuw(NumericMatrix ds_, CharacterVector penalty_, double lambda1, NumericVector lambda0s_, double theta, double aa, double bb, double eps, int max_iter, bool init_random, bool silent, bool warm_start);
RcppExport SEXP _recombine_SHC_SSL_getuw(SEXP ds_SEXP, SEXP penalty_SEXP, SEXP lambda1SEXP, SEXP lambda0s_SEXP, SEXP thetaSEXP, SEXP aaSEXP, SEXP bbSEXP, SEXP epsSEXP, SEXP max_iterSEXP, SEXP init_randomSEXP, SEXP silentSEXP, SEXP warm_startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ds_(ds_SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type penalty_(penalty_SEXP);
    Rcpp::traits::input_parameter< double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda0s_(lambda0s_SEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type aa(aaSEXP);
    Rcpp::traits::input_parameter< double >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type init_random(init_randomSEXP);
    Rcpp::traits::input_parameter< bool >::type silent(silentSEXP);
    Rcpp::traits::input_parameter< bool >::type warm_start(warm_startSEXP);
    rcpp_result_gen = Rcpp::wrap(SHC_SSL_getuw(ds_, penalty_, lambda1, lambda0s_, theta, aa, bb, eps, max_iter, init_random, silent, warm_start));
    return rcpp_result_gen;
END_RCPP
}
// SHC_get_crit
List SHC_get_crit(NumericMatrix ds_, NumericVector w_);
RcppExport SEXP _recombine_SHC_get_crit(SEXP ds_SEXP, SEXP w_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ds_(ds_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_(w_SEXP);
    rcpp_result_gen = Rcpp::wrap(SHC_get_crit(ds_, w_));
    return rcpp_result_gen;
END_RCPP
}
// SHC_get_u
List SHC_get_u(NumericMatrix ds_, NumericVector w_);
RcppExport SEXP _recombine_SHC_get_u(SEXP ds_SEXP, SEXP w_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ds_(ds_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_(w_SEXP);
    rcpp_result_gen = Rcpp::wrap(SHC_get_u(ds_, w_));
    return rcpp_result_gen;
END_RCPP
}
// multfun
IntegerMatrix multfun(IntegerMatrix x_);
RcppExport SEXP _recombine_multfun(SEXP x_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x_(x_SEXP);
    rcpp_result_gen = Rcpp::wrap(multfun(x_));
    return rcpp_result_gen;
END_RCPP
}
// distfun
NumericMatrix distfun(NumericMatrix x_);
RcppExport SEXP _recombine_distfun(SEXP x_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x_(x_SEXP);
    rcpp_result_gen = Rcpp::wrap(distfun(x_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_recombine_get_outliers_from_dist", (DL_FUNC) &_recombine_get_outliers_from_dist, 5},
    {"_recombine_get_outliers_from_ds", (DL_FUNC) &_recombine_get_outliers_from_ds, 5},
    {"_recombine_SHC_lasso_getuw", (DL_FUNC) &_recombine_SHC_lasso_getuw, 6},
    {"_recombine_SHC_lasso_lagrange_getuw", (DL_FUNC) &_recombine_SHC_lasso_lagrange_getuw, 6},
    {"_recombine_SHC_FL_getuw", (DL_FUNC) &_recombine_SHC_FL_getuw, 6},
    {"_recombine_SHC_SSL_getuw", (DL_FUNC) &_recombine_SHC_SSL_getuw, 12},
    {"_recombine_SHC_get_crit", (DL_FUNC) &_recombine_SHC_get_crit, 2},
    {"_recombine_SHC_get_u", (DL_FUNC) &_recombine_SHC_get_u, 2},
    {"_recombine_multfun", (DL_FUNC) &_recombine_multfun, 1},
    {"_recombine_distfun", (DL_FUNC) &_recombine_distfun, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_recombine(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}