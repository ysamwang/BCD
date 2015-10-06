#include "sem_el_fit.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double sem_el_fit_obj(SEXP b_weights_r, SEXP y_r, SEXP omega_r, SEXP b_r , SEXP dual_r, double tol, int max_iter) {

    //initialize graph object
  el_sem graph = el_sem(b_weights_r, y_r, omega_r, b_r, dual_r);
  return graph.update_dual(tol, max_iter);
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List sem_el_fit_weights(SEXP b_weights_r, SEXP y_r, SEXP omega_r, SEXP b_r , SEXP dual_r, double tol, int max_iter) {

    //initialize graph object
  el_sem graph = el_sem(b_weights_r, y_r, omega_r, b_r, dual_r);
  double objective = graph.update_dual(tol, max_iter);
  return Rcpp::List::create(Rcpp::Named("d", graph.get_d()),
                            Rcpp::Named("iter", graph.get_iter()),
                            Rcpp::Named("criteria", graph.get_conv_crit()),
                            Rcpp::Named("objective", objective));
}


#include "sem_el_fit.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double sem_el_naive_fit_obj(SEXP weights_r, SEXP y_r, SEXP omega_r, SEXP b_r , SEXP dual_r, double tol, int max_iter) {

    //initialize graph object
  el_sem_naive graph = el_sem_naive(weights_r, y_r, omega_r, b_r, dual_r);
  return graph.update_dual(tol, max_iter);
}

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List sem_el_naive_fit_weights(SEXP weights_r, SEXP y_r, SEXP omega_r, SEXP b_r , SEXP dual_r, double tol, int max_iter) {

    //initialize graph object
  el_sem_naive graph = el_sem_naive(weights_r, y_r, omega_r, b_r, dual_r);
  double objective = graph.update_dual(tol, max_iter);
  return Rcpp::List::create(Rcpp::Named("d", graph.get_d()),
                            Rcpp::Named("iter", graph.get_iter()),
                            Rcpp::Named("criteria", graph.get_conv_crit()),
                            Rcpp::Named("objective", objective));
}
