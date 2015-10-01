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
    el_sem(SEXP y_r, SEXP b_r, SEXP omega_r, SEXP b_weights_r, SEXP d_r,SEXP lambda_r, SEXP gamma_r, int v_r);
    void update_lambda_gamma(double tol, int max_iter);
    double get_objective();
    vec get_d();

protected:
    mat y_; // p by n data matrix
    mat omega_; // matrix of 1's and 0's representing bi-directed edge structure
    mat b_; // matrix of 1's and 0's representing directed edge structure
    mat b_weights_; // matrix of directed egdge weights
    vec d_; //vector of denominators of p_n (ie p_n = 1 / d_(n))
    double alpha_;
    vec lambda_; // dual variable for mean restrictions
    vec gamma_; //dual variable for covariance restrictions

    mat residuals_;
    uvec gamma_indices_;
    mat cross_terms_;

    int v_;
    int n_;

    void set_d_star();
    void set_gradient_hessian(vec &grad, mat &hessian, double b);
    int backtracking(vec update);

private:
};
#endif // LINSEM_H

