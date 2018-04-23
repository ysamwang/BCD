// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bcdC
Rcpp::List bcdC(SEXP Br, SEXP Omegar, SEXP BInitr, SEXP OmegaInitr, SEXP Yr, int maxIter, int sigConv, double maxKap, double tol, double omegaInitScale);
RcppExport SEXP _BCD_bcdC(SEXP BrSEXP, SEXP OmegarSEXP, SEXP BInitrSEXP, SEXP OmegaInitrSEXP, SEXP YrSEXP, SEXP maxIterSEXP, SEXP sigConvSEXP, SEXP maxKapSEXP, SEXP tolSEXP, SEXP omegaInitScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Br(BrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Omegar(OmegarSEXP);
    Rcpp::traits::input_parameter< SEXP >::type BInitr(BInitrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type OmegaInitr(OmegaInitrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Yr(YrSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< int >::type sigConv(sigConvSEXP);
    Rcpp::traits::input_parameter< double >::type maxKap(maxKapSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type omegaInitScale(omegaInitScaleSEXP);
    rcpp_result_gen = Rcpp::wrap(bcdC(Br, Omegar, BInitr, OmegaInitr, Yr, maxIter, sigConv, maxKap, tol, omegaInitScale));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BCD_bcdC", (DL_FUNC) &_BCD_bcdC, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_BCD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
