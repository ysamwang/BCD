#ifndef EL_SEM_TIE_H
#define EL_SEM_TIE_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "elSEM.h"

using namespace Rcpp ;
using namespace arma;


class el_sem_tie: public el_sem
{
public:
    //constructor
    el_sem_tie(SEXP b_weights_r, SEXP y_r, SEXP counts_r, SEXP omega_r,SEXP b_r , SEXP dual_r, int meanEst);

    double update_dual(double tol, int max_iter);
    vec get_d();

protected:
    vec counts_;

    void set_gradient_hessian(vec &grad, mat &hessian);
    int backtracking(vec update);



private:
};
#endif // EL_SEM_TIE_H
