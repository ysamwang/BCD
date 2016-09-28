#ifndef EL_SEM__NAIVE_H
#define EL_SEM__NAIVE_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "el_sem_settings.h"

using namespace Rcpp ;
using namespace arma;



class el_sem_naive
{
public:
    // Initialize sem object:
    el_sem_naive(SEXP weights_r, SEXP y_r, SEXP b_r, SEXP omega_r, bool mean_est_r);

    double update_dual(double tol, int max_iter);
    vec get_d();
    vec getGradient();

protected:
    vec d_; //vector of denominators of p_n (ie p_n = 1 / d_(n))
    vec dual_; // dual variables
    mat y_;

    bool mean_est_;
    uvec b_spots_;
    uvec omega_spots_;

    mat constraints_;

    int v_;
    int n_;

    void set_gradient_hessian(vec &grad, mat &hessian);
    int backtracking(vec update);


private:
};
#endif

