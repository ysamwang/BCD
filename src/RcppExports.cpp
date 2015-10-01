// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bcdC
Rcpp::List bcdC(SEXP Br, SEXP Omegar, SEXP BInitr, SEXP OmegaInitr, SEXP Yr, int maxIter, int sigConv, double maxKap, double tol, double omegaInitScale);
RcppExport SEXP BCD_bcdC(SEXP BrSEXP, SEXP OmegarSEXP, SEXP BInitrSEXP, SEXP OmegaInitrSEXP, SEXP YrSEXP, SEXP maxIterSEXP, SEXP sigConvSEXP, SEXP maxKapSEXP, SEXP tolSEXP, SEXP omegaInitScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
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
    __result = Rcpp::wrap(bcdC(Br, Omegar, BInitr, OmegaInitr, Yr, maxIter, sigConv, maxKap, tol, omegaInitScale));
    return __result;
END_RCPP
}
// sem_el_fitC
SEXP sem_el_fitC(SEXP y_r, SEXP b_r, SEXP omega_r, SEXP b_weights_r, SEXP d_r, SEXP lambda_r, SEXP gamma_r, int v, double tol, int max_iter);
RcppExport SEXP BCD_sem_el_fitC(SEXP y_rSEXP, SEXP b_rSEXP, SEXP omega_rSEXP, SEXP b_weights_rSEXP, SEXP d_rSEXP, SEXP lambda_rSEXP, SEXP gamma_rSEXP, SEXP vSEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type y_r(y_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b_r(b_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type omega_r(omega_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b_weights_r(b_weights_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type d_r(d_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lambda_r(lambda_rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type gamma_r(gamma_rSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    __result = Rcpp::wrap(sem_el_fitC(y_r, b_r, omega_r, b_weights_r, d_r, lambda_r, gamma_r, v, tol, max_iter));
    return __result;
END_RCPP
}
