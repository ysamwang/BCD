#include "elSEM_fixed.h"

el_sem_fixed::el_sem_fixed(SEXP b_weights_r, SEXP y_r, SEXP omega_r,SEXP b_r , SEXP dual_r, int meanEst,
                            double b_fixed, int row_ind, int col_ind) : el_sem::el_sem()
{
    //matrix which holds the observed data
    mat y = as<arma::mat>(y_r);
    n_ = y.n_cols; //number of observations
    v_ = y.n_rows; //number of variables
    counter_ = 0; //number of iterations
    conv_crit_ = 1.0; //convergence critiera (norm of gradient)

    //find appropriate spots to put in b_weights_r
    uvec b_spots = find(as<arma::mat>(b_r)); //non-structural zeros in B
    vec means = vec(v_, fill::zeros); //non-zero means
    mat b_weights = mat(v_, v_, fill::zeros); // matrix with edge weights

    if(meanEst){
        means = as<arma::vec>(b_weights_r).head(v_);
        b_weights.elem(b_spots) = as<arma::vec>(b_weights_r).tail(as<arma::vec>(b_weights_r).n_elem - v_);

    } else {
        b_weights.elem(b_spots) = as<arma::vec>(b_weights_r);
    }

    // Insert fixed element into B
    b_weights(row_ind - 1, col_ind - 1) = b_fixed;

    dual_ = as<arma::vec>(dual_r); //dual variables

    gamma_indices_ = arma::find(trimatu(as<arma::mat>(omega_r) == 0 )); //structural 0's in Omega

    constraints_ = mat(v_ + gamma_indices_.n_elem, n_);  //constraints containing mean and covariance restrictions
    constraints_.rows(0, v_ - 1) = (eye(v_, v_) - b_weights) * y; //mean restrictions

    if(meanEst) {
        constraints_.rows(0, v_ - 1).each_col() -= means;
    }


    //covariance restrictions
    int i,j,k;
    for(k = 0; k < gamma_indices_.n_elem; k++) {
        j = (int) gamma_indices_(k) / v_;
        i = (int) gamma_indices_(k) % v_;
        constraints_.row(k + v_) = constraints_.row(i) % constraints_.row(j);
    }
    d_ = (constraints_.t() * dual_) + 1.0;
}
