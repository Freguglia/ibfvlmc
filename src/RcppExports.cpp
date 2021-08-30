// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ibf
List ibf(List z_test, IntegerVector z_train, IntegerVector renewal, LogicalMatrix allowedMatrix, double alpha, double logprior_penalty, unsigned int Hmax, unsigned int alphlen, unsigned int burnin, unsigned int nsamples);
RcppExport SEXP _ibfvlmc_ibf(SEXP z_testSEXP, SEXP z_trainSEXP, SEXP renewalSEXP, SEXP allowedMatrixSEXP, SEXP alphaSEXP, SEXP logprior_penaltySEXP, SEXP HmaxSEXP, SEXP alphlenSEXP, SEXP burninSEXP, SEXP nsamplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type z_test(z_testSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z_train(z_trainSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type renewal(renewalSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type allowedMatrix(allowedMatrixSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type logprior_penalty(logprior_penaltySEXP);
    Rcpp::traits::input_parameter< unsigned int >::type Hmax(HmaxSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type alphlen(alphlenSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nsamples(nsamplesSEXP);
    rcpp_result_gen = Rcpp::wrap(ibf(z_test, z_train, renewal, allowedMatrix, alpha, logprior_penalty, Hmax, alphlen, burnin, nsamples));
    return rcpp_result_gen;
END_RCPP
}
// ibf_comp
List ibf_comp(List z_test, IntegerVector z_train, IntegerVector renewal, LogicalMatrix allowedMatrix, double alpha, double logprior_penalty, unsigned int Hmax, unsigned int alphlen, unsigned int burnin, unsigned int nsamples);
RcppExport SEXP _ibfvlmc_ibf_comp(SEXP z_testSEXP, SEXP z_trainSEXP, SEXP renewalSEXP, SEXP allowedMatrixSEXP, SEXP alphaSEXP, SEXP logprior_penaltySEXP, SEXP HmaxSEXP, SEXP alphlenSEXP, SEXP burninSEXP, SEXP nsamplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type z_test(z_testSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z_train(z_trainSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type renewal(renewalSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type allowedMatrix(allowedMatrixSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type logprior_penalty(logprior_penaltySEXP);
    Rcpp::traits::input_parameter< unsigned int >::type Hmax(HmaxSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type alphlen(alphlenSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nsamples(nsamplesSEXP);
    rcpp_result_gen = Rcpp::wrap(ibf_comp(z_test, z_train, renewal, allowedMatrix, alpha, logprior_penalty, Hmax, alphlen, burnin, nsamples));
    return rcpp_result_gen;
END_RCPP
}
// rvlmc_cpp
IntegerVector rvlmc_cpp(unsigned int n, List context_list, List probs);
RcppExport SEXP _ibfvlmc_rvlmc_cpp(SEXP nSEXP, SEXP context_listSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< List >::type context_list(context_listSEXP);
    Rcpp::traits::input_parameter< List >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(rvlmc_cpp(n, context_list, probs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ibfvlmc_ibf", (DL_FUNC) &_ibfvlmc_ibf, 10},
    {"_ibfvlmc_ibf_comp", (DL_FUNC) &_ibfvlmc_ibf_comp, 10},
    {"_ibfvlmc_rvlmc_cpp", (DL_FUNC) &_ibfvlmc_rvlmc_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ibfvlmc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
