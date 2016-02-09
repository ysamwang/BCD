#include "elSEM_fixed.h"

el_sem_fixed::el_sem_fixed(SEXP b_weights_r, SEXP y_r,SEXP b_r, SEXP moment_2_restrictions_r, SEXP moment_3_restrictions_r, int meanEst,
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

 // Initialize the matrix of constraints and the vector of dual variables (start at 0)
    mat moment_2_restrictions = as<arma::mat>(moment_2_restrictions_r);
    mat moment_3_restrictions;

    if(!Rf_isNull(moment_3_restrictions_r)) {
        moment_3_restrictions = as<arma::mat>(moment_3_restrictions_r);
        constraints_ = mat(v_ + moment_2_restrictions.n_rows + moment_3_restrictions.n_rows, n_);  //constraints containing mean and covariance  and 3rd moment restrictions

    } else {
        constraints_ = mat(v_ + moment_2_restrictions.n_rows, n_);  //constraints containing mean and covariance restrictions
    }

    dual_ = vec(constraints_.n_rows, fill::zeros);


    // Filling in the constraint matrix

    // Mean restrictions
    constraints_.rows(0, v_ - 1) = (eye(v_, v_) - b_weights) * y;
    // if there is a non-zero mean vector
    if(meanEst) {
        constraints_.rows(0, v_ - 1).each_col() -= means;
    }



    // covariance restrictions
    int k;
    for(k = 0; k < moment_2_restrictions.n_rows; k++) {
            // 2nd moment constraint
            constraints_.row(k + v_) = constraints_.row( moment_2_restrictions(k, 0) ) % constraints_.row( moment_2_restrictions(k, 1) );
        }

    //3rd moment restrictions if needed
    if(!Rf_isNull(moment_3_restrictions_r)) {
        for(k = 0; k < moment_3_restrictions.n_rows; k++){
            constraints_.row(k + v_ + moment_2_restrictions.n_rows) = constraints_.row( moment_3_restrictions(k, 0) ) % constraints_.row( moment_3_restrictions(k, 1) ) % constraints_.row( moment_3_restrictions(k, 2) );
        }
    }

    d_ = (constraints_.t() * dual_) + 1.0;
}
