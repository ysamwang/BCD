#ifndef BCD_H
#define BCD_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "linSEM.h"

using namespace Rcpp ;
using namespace arma;

Rcpp::List bcdC(Rcpp::List model_r, int maxIter, int sigConv, double maxKap, double tol);


#endif // LINSEM_H

