#include "sempl.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
NumericVector sempl_input(SEXP b_weights_r, SEXP y_r, SEXP b_r, bool mean_est_r, SEXP covar_restrict,
                          double tol, int max_iter)
{
    // squash warnings about singular hessian
    if(1) {
        std::ostream nullstream(0);
        set_stream_err2(nullstream);
    }
    //initialize graph object
    el_sem graph = el_sem(b_weights_r, y_r, b_r, mean_est_r, covar_restrict);

    NumericVector ret = NumericVector::create(-graph.update_dual(tol, max_iter));
    ret.attr("gradient") = -graph.getGradient();
    return ret;
}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List sempl_input_weights(SEXP b_weights_r, SEXP y_r, SEXP b_r, bool mean_est_r, SEXP covar_restrict,
                               double tol, int max_iter)
{
    // squash warnings about singular hessian
    if(1) {
        std::ostream nullstream(0);
        set_stream_err2(nullstream);
    }
    //initialize graph object
    el_sem graph = el_sem(b_weights_r, y_r, b_r, mean_est_r, covar_restrict);

    double ret = -graph.update_dual(tol, max_iter);

    return Rcpp::List::create(Rcpp::Named("d", graph.get_d()),
                              Rcpp::Named("lr", ret));
}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
NumericVector sempl_input_naive(SEXP weights_r, SEXP y_r, SEXP b_r, SEXP omega_r, bool mean_est_r,
                                double tol, int max_iter)
{
    // squash warnings about singular hessian
    if(1) {
        std::ostream nullstream(0);
        set_stream_err2(nullstream);
    }
    //initialize graph object
    el_sem_naive graph = el_sem_naive(weights_r, y_r, b_r, omega_r, mean_est_r);

    NumericVector ret = NumericVector::create(-graph.update_dual(tol, max_iter));
    ret.attr("gradient") = -graph.getGradient();
    return ret;
}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List sempl_input_naive_weights(SEXP weights_r, SEXP y_r, SEXP b_r, SEXP omega_r, bool mean_est_r,
                                     double tol, int max_iter)
{
    // squash warnings about singular hessian
    if(1) {
        std::ostream nullstream(0);
        set_stream_err2(nullstream);
    }
    //initialize graph object
    el_sem_naive graph = el_sem_naive(weights_r, y_r, b_r, omega_r, mean_est_r);

    NumericVector ret = NumericVector::create(-graph.update_dual(tol, max_iter));
    return Rcpp::List::create(Rcpp::Named("d", graph.get_d()),
                              Rcpp::Named("lr", ret));
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
NumericVector sempl_input_grad(SEXP b_weights_r, SEXP y_r, SEXP b_r, bool mean_est_r, SEXP covar_restrict,
                          double tol, int max_iter)
{
    // squash warnings about singular hessian
    if(1) {
        std::ostream nullstream(0);
        set_stream_err2(nullstream);
    }
    //initialize graph object
    el_sem graph = el_sem(b_weights_r, y_r, b_r, mean_est_r, covar_restrict);
    graph.update_dual(tol, max_iter);
    vec grad = -graph.getGradient();
    NumericVector ret = Rcpp::as<NumericVector>(wrap(grad));
    return ret;
}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double sempl_input_reg(SEXP b_weights_r, SEXP y_r, SEXP b_r, bool mean_est_r, SEXP covar_restrict,
                          double tol, int max_iter)
{
    // squash warnings about singular hessian
    if(1) {
        std::ostream nullstream(0);
        set_stream_err2(nullstream);
    }
    //initialize graph object

    el_sem graph = el_sem(b_weights_r, y_r, b_r, mean_est_r, covar_restrict);

    double ret = -graph.update_dual(tol, max_iter);
    return ret;
}
