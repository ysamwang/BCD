#include "sem_el_fit.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double sem_el_fitC(SEXP y_r, SEXP omega_r, SEXP b_weights_r, SEXP dual_r
                   , int v, double tol, int max_iter) {

    //initialize graph object
  el_sem graph = el_sem(y_r, omega_r, b_weights_r, d_r, dual_r, v);
  return graph.update_dual(tol, max_iter);
}
