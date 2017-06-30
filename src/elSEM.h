#ifndef EL_SEM_H
#define EL_SEM_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>
#include "el_sem_settings.h"

using namespace Rcpp ;
using namespace arma;

class el_sem
{
public:
    // Initialize sem object:
    el_sem(SEXP b_weights_r, SEXP y_r, SEXP b_r, bool mean_est_r, SEXP covar_restrict);


    double update_dual(double tol, int max_iter);
    vec get_d();
    vec getGradient();

protected:
    mat y_;
    vec d_; //vector of denominators of p_n (ie p_n = 1 / d_(n))
    vec dual_; // dual variables

    mat constraints_;
    mat restrict_;

    int v_;
    int n_;
    int counter_;
    double conv_crit_;
    bool mean_est_;

    uvec b_spots_;

    void set_gradient_hessian(vec &grad, mat &hessian);
    int backtracking(vec update);


private:
};
#endif // ELSEM_H

