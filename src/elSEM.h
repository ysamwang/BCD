#ifndef EL_SEM_H
#define EL_SEM_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp ;
using namespace arma;

class el_sem
{
public:
    // Initialize sem object:
    el_sem(SEXP b_weights_r, SEXP y_r, SEXP omega_r,SEXP b_r , SEXP dual_r, int meanEst);
    el_sem();

    double update_dual(double tol, int max_iter);
    vec get_d();

    double get_conv_crit();
    int get_iter();

protected:
    vec d_; //vector of denominators of p_n (ie p_n = 1 / d_(n))
    vec dual_; // dual variables

    uvec gamma_indices_;
    mat constraints_;

    int v_;
    int n_;
    int counter_;
    double conv_crit_;

    void set_gradient_hessian(vec &grad, mat &hessian);
    int backtracking(vec update);

private:
};
#endif // ELSEM_H

