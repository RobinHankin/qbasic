// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pochhammer
NumericVector pochhammer(const NumericVector a_real, const NumericVector a_imag, const NumericVector q_real, const NumericVector q_imag, const NumericVector n, const NumericVector maxit);
RcppExport SEXP _queueR_pochhammer(SEXP a_realSEXP, SEXP a_imagSEXP, SEXP q_realSEXP, SEXP q_imagSEXP, SEXP nSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type a_real(a_realSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type a_imag(a_imagSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type q_real(q_realSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type q_imag(q_imagSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(pochhammer(a_real, a_imag, q_real, q_imag, n, maxit));
    return rcpp_result_gen;
END_RCPP
}
// qexp_C
NumericVector qexp_C(const NumericVector z_real, const NumericVector z_imag, const NumericVector q_real, const NumericVector q_imag, const NumericVector maxit);
RcppExport SEXP _queueR_qexp_C(SEXP z_realSEXP, SEXP z_imagSEXP, SEXP q_realSEXP, SEXP q_imagSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type z_real(z_realSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type z_imag(z_imagSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type q_real(q_realSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type q_imag(q_imagSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(qexp_C(z_real, z_imag, q_real, q_imag, maxit));
    return rcpp_result_gen;
END_RCPP
}
// BasicHypergeometric
NumericVector BasicHypergeometric(const NumericVector z_real, const NumericVector z_imag, const NumericVector q_real, const NumericVector q_imag, const NumericVector a_real, const NumericVector a_imag, const NumericVector b_real, const NumericVector b_imag, const NumericVector maxit);
RcppExport SEXP _queueR_BasicHypergeometric(SEXP z_realSEXP, SEXP z_imagSEXP, SEXP q_realSEXP, SEXP q_imagSEXP, SEXP a_realSEXP, SEXP a_imagSEXP, SEXP b_realSEXP, SEXP b_imagSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type z_real(z_realSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type z_imag(z_imagSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type q_real(q_realSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type q_imag(q_imagSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type a_real(a_realSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type a_imag(a_imagSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type b_real(b_realSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type b_imag(b_imagSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(BasicHypergeometric(z_real, z_imag, q_real, q_imag, a_real, a_imag, b_real, b_imag, maxit));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_queueR_pochhammer", (DL_FUNC) &_queueR_pochhammer, 6},
    {"_queueR_qexp_C", (DL_FUNC) &_queueR_qexp_C, 5},
    {"_queueR_BasicHypergeometric", (DL_FUNC) &_queueR_BasicHypergeometric, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_queueR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
