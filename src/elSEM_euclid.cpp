#include "elSEM_euclid.h"

double el_sem_euclid_obj(SEXP b_weights_r, SEXP y_r, SEXP omega_r, SEXP b_r)
{
    //matrix which holds the observed data
    mat y = as<arma::mat>(y_r);
    int n = y.n_cols; //number of observations
    int v = y.n_rows; //number of variables

    //find appropriate spots to put in b_weights_r
    uvec b_spots = find(as<arma::mat>(b_r)); //non-structural zeros in B
    mat b_weights(v, v, fill::zeros); // matrix with edge weights
    b_weights.elem(b_spots) = as<arma::vec>(b_weights_r);

    uvec gamma_indices = arma::find(trimatu(as<arma::mat>(omega_r) == 0 )); //structural 0's in Omega

    mat constraints(v + gamma_indices.n_elem, n); //constraints containing mean and covariance restrictions
    constraints.rows(0, v - 1) = (eye(v, v) - b_weights) * y ; //mean restrictions

    //covariance restrictions
    int i,j,k;
    for(k = 0; k < gamma_indices.n_elem; k++) {
        j = (int) gamma_indices(k) / v;
        i = (int) gamma_indices(k) % v;
        constraints.row(k + v) = constraints.row(i) % constraints.row(j);
    }

    vec avg_constraints = mean(constraints, 1);

    mat S(v + gamma_indices.n_elem, v + gamma_indices.n_elem, fill::zeros);
    for(i = 0; i < n; i++){
        S += constraints.col(i) * (avg_constraints.t() - constraints.col(i).t() );
    }
    S = S / n;

    vec dual = solve(S, avg_constraints);

    constraints.each_col() -= avg_constraints;
    vec p_star = (1.0/ n) * (1 + (dual.t() * constraints).t());


    double objective = - 1.0/2 * sum(pow(n * p_star - 1.0 ,2));

    return objective;
}

Rcpp::List el_sem_euclid_weights(SEXP b_weights_r, SEXP y_r, SEXP omega_r, SEXP b_r)
{
    //matrix which holds the observed data
    mat y = as<arma::mat>(y_r);
    int n = y.n_cols; //number of observations
    int v = y.n_rows; //number of variables

    //find appropriate spots to put in b_weights_r
    uvec b_spots = find(as<arma::mat>(b_r)); //non-structural zeros in B
    mat b_weights(v, v, fill::zeros); // matrix with edge weights
    b_weights.elem(b_spots) = as<arma::vec>(b_weights_r);

    uvec gamma_indices = arma::find(trimatu(as<arma::mat>(omega_r) == 0 )); //structural 0's in Omega

    mat constraints(v + gamma_indices.n_elem, n); //constraints containing mean and covariance restrictions
    constraints.rows(0, v - 1) = (eye(v, v) - b_weights) * y ; //mean restrictions




    //covariance restrictions
    int i,j,k;
    for(k = 0; k < gamma_indices.n_elem; k++) {
        j = (int) gamma_indices(k) / v;
        i = (int) gamma_indices(k) % v;
        constraints.row(k + v) = constraints.row(i) % constraints.row(j);
    }

    vec avg_constraints = mean(constraints, 1);
    mat S(v + gamma_indices.n_elem, v + gamma_indices.n_elem, fill::zeros);
    for(i = 0; i < n; i++){
        S += constraints.col(i) * (avg_constraints.t() - constraints.col(i).t() );
    }
    S = S / n;

    vec dual = solve(S, avg_constraints);

    constraints.each_col() -= avg_constraints;
    vec p_star = (1.0/ n) * (1 + (dual.t() * constraints).t());


    double objective = - 1.0/2 * sum(pow(n * p_star - 1.0 ,2));

    return Rcpp::List::create(Rcpp::Named("p_star", p_star),
                            Rcpp::Named("objective", objective));
}
