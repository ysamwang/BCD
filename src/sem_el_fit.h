#ifndef SEM_EL_FIT_H
#define SEM_EL_FIT_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "elSEM.h"

using namespace Rcpp ;
using namespace arma;

double sem_el_fitC(SEXP y_r, SEXP b_r, SEXP omega_r, SEXP b_weights_r, SEXP omega_weights_r,
                    SEXP lambda_r, SEXP gamma_r, double tol, int max_iter);


#endif // LINSEM_H

