#include "sem_el_fit.h"


// Main Functions

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double sem_el_fit_obj(SEXP b_weights_r, SEXP y_r, SEXP b_r,  SEXP moment_2_restrictions_r, SEXP moment_3_restrictions_r, SEXP moment_4_restrictions_r, double tol, int max_iter, int meanEst)
{
    // squash warnings about singular hessian
    if(1){
        std::ostream nullstream(0);
        set_stream_err2(nullstream);
    }
    //initialize graph object
    el_sem graph = el_sem(b_weights_r, y_r, b_r, moment_2_restrictions_r, moment_3_restrictions_r, moment_4_restrictions_r, meanEst);

    double ret = graph.update_dual(tol, max_iter);
    return ret;
}


double sem_el_fit_obj_nlm(SEXP b_weights_r, SEXP y_r, SEXP b_r,  SEXP moment_2_restrictions_r, SEXP moment_3_restrictions_r, SEXP moment_4_restrictions_r, double tol, int max_iter, int meanEst)
{
  // squash warnings about singular hessian
  if(1){
    std::ostream nullstream(0);
    set_stream_err2(nullstream);
  }
  //initialize graph object
  el_sem graph = el_sem(b_weights_r, y_r, b_r, moment_2_restrictions_r, moment_3_restrictions_r, moment_4_restrictions_r, meanEst);
  
  NumericVector ret = NumericVector::create(graph.update_dual(tol, max_iter));
  ret.attr("gradient") = graph.getGradient();
  return ret;
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List sem_el_fit_weights(SEXP b_weights_r, SEXP y_r, SEXP b_r,  SEXP moment_2_restrictions_r, SEXP moment_3_restrictions_r, SEXP moment_4_restrictions_r, double tol, int max_iter, int meanEst)
{
    //initialize graph object
    el_sem graph = el_sem(b_weights_r, y_r, b_r, moment_2_restrictions_r, moment_3_restrictions_r, moment_4_restrictions_r, meanEst);
    double objective = graph.update_dual(tol, max_iter);
    return Rcpp::List::create(Rcpp::Named("d", graph.get_d()),
                              Rcpp::Named("iter", graph.get_iter()),
                              Rcpp::Named("dual", graph.get_dual()),
                              Rcpp::Named("criteria", graph.get_conv_crit()),
                              Rcpp::Named("objective", objective));
}

// Naive Functions


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double sem_el_naive_fit_obj(SEXP weights_r, SEXP y_r, SEXP omega_r, SEXP b_r , SEXP dual_r, double tol, int max_iter, int meanEst)
{

    //initialize graph object
    el_sem_naive graph = el_sem_naive(weights_r, y_r, omega_r, b_r, dual_r, meanEst);
    return graph.update_dual(tol, max_iter);
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List sem_el_naive_fit_weights(SEXP weights_r, SEXP y_r, SEXP omega_r, SEXP b_r , SEXP dual_r, double tol, int max_iter, int meanEst)
{

    //initialize graph object
    el_sem_naive graph = el_sem_naive(weights_r, y_r, omega_r, b_r, dual_r, meanEst);
    double objective = graph.update_dual(tol, max_iter);
    return Rcpp::List::create(Rcpp::Named("d", graph.get_d()),
                              Rcpp::Named("iter", graph.get_iter()),
                              Rcpp::Named("criteria", graph.get_conv_crit()),
                              Rcpp::Named("objective", objective));
}


// Fixed Functions

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double sem_el_fit_obj_one_fixed(SEXP b_weights_r, SEXP y_r, SEXP b_r, SEXP moment_2_restrictions_r, SEXP moment_3_restrictions_r, double tol, int max_iter, int meanEst,
                                double b_fixed, int row_ind, int col_ind)
{

    //initialize graph object
    el_sem_fixed graph = el_sem_fixed(b_weights_r, y_r, b_r, moment_2_restrictions_r, moment_3_restrictions_r, meanEst,
                                      b_fixed, row_ind, col_ind);

    double ret = graph.update_dual(tol, max_iter);
    return ret;
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double sem_el_euclid_fit_obj(SEXP b_weights_r, SEXP y_r, SEXP b_r, SEXP moment_2_restrictions_r,
                              SEXP moment_3_restrictions_r,  SEXP moment_4_restrictions_r, int meanEst)
{
        // squash warnings about singular hessian
    if(1){
        std::ostream nullstream(0);
        set_stream_err2(nullstream);
    }
    return el_sem_euclid_obj(b_weights_r, y_r, b_r, moment_2_restrictions_r,
                              moment_3_restrictions_r, moment_4_restrictions_r, meanEst);
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List sem_el_euclid_fit_weights(SEXP b_weights_r, SEXP y_r, SEXP b_r, SEXP moment_2_restrictions_r,
                              SEXP moment_3_restrictions_r,  SEXP moment_4_restrictions_r, int meanEst)
{
        // squash warnings about singular hessian
    if(1){
        std::ostream nullstream(0);
        set_stream_err2(nullstream);
    }
    return el_sem_euclid_weights(b_weights_r, y_r, b_r, moment_2_restrictions_r,
                              moment_3_restrictions_r, moment_4_restrictions_r, meanEst);
}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double sem_el_fit_obj_l1(SEXP b_weights_r, SEXP y_r, SEXP b_r,  SEXP moment_2_restrictions_r, SEXP moment_3_restrictions_r, SEXP moment_4_restrictions_r, double tol, int max_iter, int meanEst, double penalty)
{
    // squash warnings about singular hessian
    if(1){
        std::ostream nullstream(0);
        set_stream_err2(nullstream);
    }
    //initialize graph object
    el_sem graph = el_sem(b_weights_r, y_r, b_r, moment_2_restrictions_r, moment_3_restrictions_r, moment_4_restrictions_r, meanEst);
    double ret;
    if(meanEst) {
        int v_ = as<arma::mat>(y_r).n_rows;
        ret = graph.update_dual(tol, max_iter) + penalty * accu(abs(as<arma::vec>(b_weights_r).tail(as<arma::vec>(b_weights_r).n_elem - v_)));

    } else {
        ret = graph.update_dual(tol, max_iter) + penalty * accu(abs(as<arma::vec>(b_weights_r)));
    }

    return ret;
}
