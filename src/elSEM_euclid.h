#ifndef EL_SEM_EUCLID_H
#define EL_SEM_EUCLID_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp ;
using namespace arma;


Rcpp::List el_sem_euclid_weights(SEXP b_weights_r, SEXP y_r, SEXP omega_r, SEXP b_r);
double el_sem_euclid_obj(SEXP b_weights_r, SEXP y_r, SEXP omega_r, SEXP b_r);

#endif // ELSEM_H

