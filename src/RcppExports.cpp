// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ibf
void ibf(IntegerVector z_test, List z_train, IntegerVector renewal, unsigned int Hmax, unsigned int alphlen, unsigned int burnin);
RcppExport SEXP _ibfvlmc_ibf(SEXP z_testSEXP, SEXP z_trainSEXP, SEXP renewalSEXP, SEXP HmaxSEXP, SEXP alphlenSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type z_test(z_testSEXP);
    Rcpp::traits::input_parameter< List >::type z_train(z_trainSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type renewal(renewalSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type Hmax(HmaxSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type alphlen(alphlenSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burnin(burninSEXP);
    ibf(z_test, z_train, renewal, Hmax, alphlen, burnin);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ibfvlmc_ibf", (DL_FUNC) &_ibfvlmc_ibf, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_ibfvlmc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
