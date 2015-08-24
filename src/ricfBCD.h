#ifndef RICFBCD_H
#define RICFBCD_H

#define BOOST_DISABLE_ASSERTS true

#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp ;
using namespace arma;
Rcpp::List mStep_C(mm_model model, double elbo_T, int stepType, int maxAlphaIter, int maxThetaIter, int maxLSIter,
                              double alphaTol, double thetaTol, double aNaught, double tau,
                               int bMax, double bNaught, double bMult, int vCutoff, NumericVector holdConst, NumericVector iterReached);

#endif
