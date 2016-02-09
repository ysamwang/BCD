#include "elSEM_euclid.h"

double el_sem_euclid_obj(SEXP b_weights_r, SEXP y_r, SEXP b_r, SEXP moment_2_restrictions_r, SEXP moment_3_restrictions_r,  SEXP moment_4_restrictions_r, int meanEst)
{
    //matrix which holds the observed data
    mat y = as<arma::mat>(y_r);
    int n_ = y.n_cols; //number of observations
    int v_ = y.n_rows; //number of variables

    mat constraints_;

    //find appropriate spots to put in b_weights_r
    uvec b_spots = find(as<arma::mat>(b_r)); //non-structural zeros in B
    vec means = vec(v_, fill::zeros); //non-zero means
    mat b_weights = mat(v_, v_, fill::zeros); // matrix with edge weights

    if(meanEst) {
        means = as<arma::vec>(b_weights_r).head(v_);
        b_weights.elem(b_spots) = as<arma::vec>(b_weights_r).tail(as<arma::vec>(b_weights_r).n_elem - v_);

    } else {
        b_weights.elem(b_spots) = as<arma::vec>(b_weights_r);
    }


    // Initialize the matrix of constraints and the vector of dual variables (start at 0)
    mat moment_2_restrictions = as<arma::mat>(moment_2_restrictions_r);
    mat moment_3_restrictions;
    mat moment_4_restrictions;

    if(!Rf_isNull(moment_4_restrictions_r)) {
        moment_4_restrictions = as<arma::mat>(moment_4_restrictions_r);
        moment_3_restrictions = as<arma::mat>(moment_3_restrictions_r);
        constraints_ = mat(v_ + moment_2_restrictions.n_rows + moment_3_restrictions.n_rows + moment_4_restrictions.n_rows, n_);  //constraints containing mean and covariance  and 3rd moment restrictions
    } else if(!Rf_isNull(moment_3_restrictions_r)) {
        moment_3_restrictions = as<arma::mat>(moment_3_restrictions_r);
        constraints_ = mat(v_ + moment_2_restrictions.n_rows + moment_3_restrictions.n_rows, n_);  //constraints containing mean and covariance  and 3rd moment restrictions
    } else {
        constraints_ = mat(v_ + moment_2_restrictions.n_rows, n_);  //constraints containing mean and covariance restrictions
    }

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
            for(k = 0; k < moment_3_restrictions.n_rows; k++) {
                constraints_.row(k + v_ + moment_2_restrictions.n_rows) = constraints_.row( moment_3_restrictions(k, 0) ) %
                        constraints_.row( moment_3_restrictions(k, 1) ) % constraints_.row( moment_3_restrictions(k, 2) );
            }
        }

        if(!Rf_isNull(moment_4_restrictions_r)) {
            for(k = 0; k < moment_4_restrictions.n_rows; k++) {
                constraints_.row(k + v_ + moment_2_restrictions.n_rows + moment_3_restrictions.n_rows) = constraints_.row( moment_4_restrictions(k, 0) ) %
                        constraints_.row( moment_4_restrictions(k, 1) ) % constraints_.row( moment_4_restrictions(k, 2) ) % constraints_.row( moment_4_restrictions(k, 3) );
            }
        }


    vec avg_constraints = mean(constraints_, 1);

    int i;
    mat S(constraints_.n_rows, constraints_.n_rows, fill::zeros);
    for(i = 0; i < n_; i++){
        S += (constraints_.col(i) - avg_constraints) * (constraints_.col(i) - avg_constraints).t();
    }

    S = S / n_;

    vec dual = solve(S, avg_constraints);

    constraints_.each_col() -= avg_constraints;
    vec p_star = (1.0/ n_) * (1 - (dual.t() * constraints_).t());


    double objective = - 1.0/2.0 * sum(pow(n_ * p_star - 1.0 ,2));

    return objective;
}

Rcpp::List el_sem_euclid_weights(SEXP b_weights_r, SEXP y_r, SEXP b_r, SEXP moment_2_restrictions_r, SEXP moment_3_restrictions_r,  SEXP moment_4_restrictions_r, int meanEst)
{
    //matrix which holds the observed data
    mat y = as<arma::mat>(y_r);
    int n_ = y.n_cols; //number of observations
    int v_ = y.n_rows; //number of variables

    mat constraints_;

    //find appropriate spots to put in b_weights_r
    uvec b_spots = find(as<arma::mat>(b_r)); //non-structural zeros in B
    vec means = vec(v_, fill::zeros); //non-zero means
    mat b_weights = mat(v_, v_, fill::zeros); // matrix with edge weights

    if(meanEst) {
        means = as<arma::vec>(b_weights_r).head(v_);
        b_weights.elem(b_spots) = as<arma::vec>(b_weights_r).tail(as<arma::vec>(b_weights_r).n_elem - v_);

    } else {
        b_weights.elem(b_spots) = as<arma::vec>(b_weights_r);
    }


    // Initialize the matrix of constraints and the vector of dual variables (start at 0)
    mat moment_2_restrictions = as<arma::mat>(moment_2_restrictions_r);
    mat moment_3_restrictions;
    mat moment_4_restrictions;

    if(!Rf_isNull(moment_4_restrictions_r)) {
        moment_4_restrictions = as<arma::mat>(moment_4_restrictions_r);
        moment_3_restrictions = as<arma::mat>(moment_3_restrictions_r);
        constraints_ = mat(v_ + moment_2_restrictions.n_rows + moment_3_restrictions.n_rows + moment_4_restrictions.n_rows, n_);  //constraints containing mean and covariance  and 3rd moment restrictions
    } else if(!Rf_isNull(moment_3_restrictions_r)) {
        moment_3_restrictions = as<arma::mat>(moment_3_restrictions_r);
        constraints_ = mat(v_ + moment_2_restrictions.n_rows + moment_3_restrictions.n_rows, n_);  //constraints containing mean and covariance  and 3rd moment restrictions
    } else {
        constraints_ = mat(v_ + moment_2_restrictions.n_rows, n_);  //constraints containing mean and covariance restrictions
    }

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
            for(k = 0; k < moment_3_restrictions.n_rows; k++) {
                constraints_.row(k + v_ + moment_2_restrictions.n_rows) = constraints_.row( moment_3_restrictions(k, 0) ) %
                        constraints_.row( moment_3_restrictions(k, 1) ) % constraints_.row( moment_3_restrictions(k, 2) );
            }
        }

        if(!Rf_isNull(moment_4_restrictions_r)) {
            for(k = 0; k < moment_4_restrictions.n_rows; k++) {
                constraints_.row(k + v_ + moment_2_restrictions.n_rows + moment_3_restrictions.n_rows) = constraints_.row( moment_4_restrictions(k, 0) ) %
                        constraints_.row( moment_4_restrictions(k, 1) ) % constraints_.row( moment_4_restrictions(k, 2) ) % constraints_.row( moment_4_restrictions(k, 3) );
            }
        }


    vec avg_constraints = mean(constraints_, 1);

    int i;
    mat S(constraints_.n_rows, constraints_.n_rows, fill::zeros);
    for(i = 0; i < n_; i++){
        S += (constraints_.col(i) - avg_constraints) * (constraints_.col(i) - avg_constraints).t();
    }

    S = S / n_;

    vec dual = solve(S, avg_constraints);

    constraints_.each_col() -= avg_constraints;
    vec p_star = (1.0/ n_) * (1 - (dual.t() * constraints_).t());


    double objective = - 1.0/2 * sum(pow(n_ * p_star - 1.0 ,2));

    return Rcpp::List::create(Rcpp::Named("p_star", p_star),
                            Rcpp::Named("objective", objective));
}
