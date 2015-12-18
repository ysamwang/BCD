#ifndef EL_SEM_FIXED_H
#define EL_SEM_FIXED_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "elSEM.h"

using namespace Rcpp ;
using namespace arma;


class el_sem_fixed: public el_sem
{
public:
    el_sem_fixed(SEXP b_weights_r, SEXP y_r, SEXP omega_r,SEXP b_r , SEXP dual_r, int meanEst, double b_fixed, int row_ind, int col_ind);
    el_sem_fixed(SEXP b_weights_r, SEXP y_r, SEXP omega_r,SEXP b_r , SEXP dual_r, int meanEst);

};
#endif // EL_SEM_FIXED_H
