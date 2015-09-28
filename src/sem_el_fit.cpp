#include "sem_el_fit.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
double sem_el_fitC(SEXP y_r, SEXP b_r, SEXP omega_r, SEXP b_weights_r, SEXP omega_weights_r, SEXP d_r, SEXP lambda_r,
                    SEXP gamma_r, int v, double tol, int max_iter) {

    //initialize graph object
  el_sem graph = el_sem(y_r, b_r, omega_r, b_weights_r, omega_weights_r, d_r, lambda_r, gamma_r, v);

  graph.update_lambda_gamma(tol, max_iter);
  return graph.get_objective();
}
