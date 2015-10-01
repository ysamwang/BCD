#include "elSEM.h"

el_sem::el_sem(SEXP y_r, SEXP b_r, SEXP omega_r, SEXP b_weights_r, SEXP d_r, SEXP lambda_r, SEXP gamma_r, int v_r)
{
    y_ = as<arma::mat>(y_r);
    n_ = y_.n_cols;
    b_ = as<arma::mat>(b_r);
    omega_ = as<arma::mat>(omega_r);
    b_weights_ = as<arma::mat>(b_weights_r);
    d_ = as<arma::vec>(d_r);
    lambda_ = as<arma::vec>(lambda_r);
    gamma_ = as<arma::vec>(gamma_r);
    v_ = v_r;

    gamma_indices_ = arma::find(trimatu(omega_ == 0 ));
    Rcout <<trimatu(omega_ == 0 ) <<endl;
    residuals_ = (eye(v_, v_) - b_weights_) * y_ ;
    cross_terms_.zeros(gamma_indices_.n_elem, n_);

    int i,j,k;
    for(k = 0; k < gamma_indices_.n_elem; k++) {
        j = (int) gamma_indices_(k) / v_;
        i = (int) gamma_indices_(k) % v_;
        cross_terms_.row(k) = residuals_.row(i) % residuals_.row(j);
    }
}

void el_sem::update_lambda_gamma(double tol, int max_iter)
{
    double conv_crit = 1.0;
    int counter = 0;
    double old_obj, new_obj;
    new_obj = get_objective();
    double b = 100000.0;
    double backtracking_scaling_const = .8;
    vec ret(max_iter);

    //pre-allocate memory for gradient, hessian and update
    vec grad(v_ + gamma_indices_.n_elem );
    mat hessian(v_ + gamma_indices_.n_elem , v_ + gamma_indices_.n_elem );
    vec update(v_ + gamma_indices_.n_elem );
    grad.zeros();
    hessian.zeros();
    update.zeros();

    while(conv_crit > tol && counter < max_iter) {


        // build in backtracking if necessary
        set_gradient_hessian(grad, hessian, b);

//        Rcout << "Hessian: " <<endl << hessian <<endl;
//        Rcout << "Grad: " <<endl <<grad <<endl;
//        Rcout << cond(hessian)<< std::endl;


        // back tracking line search

        update =  join_cols(lambda_, gamma_) - solve(hessian, grad);
        int counterB = 0;
        while(!backtracking(update)){
            update *= backtracking_scaling_const;
            counterB ++;
        }

        lambda_ = update.rows(0, v_ - 1);
        gamma_ = update.rows(v_, v_ + gamma_indices_.n_elem - 1);
        set_d_star();
//        conv_crit = fabs((new_obj - old_obj)/old_obj);
        conv_crit = norm( grad ,2);
        counter++;
        Rcout << "Iter: " << counter << " Backtracks: " <<counterB <<endl;
        Rcout  << "Sum: "  << sum(1.0 / d_) <<std::endl;
        Rcout << d_ <<std::endl;
        Rcout  << "Conv Criteria: " << conv_crit <<std::endl;
        Rcout << "Objective: " << get_objective() <<std::endl;

    }
}

void el_sem::set_d_star()
{
    d_.ones();
    // term from
    d_ += (lambda_.t() * residuals_).t();
    d_ += (gamma_.t() * cross_terms_).t();
    d_ *= 1.0 / n_;
 }


int el_sem::backtracking(vec update)
{
    vec d = d_;
    d.ones();
    // term from
    d += (update.subvec(0, v_ - 1).t() * residuals_).t();
    d += (update.subvec(v_, v_ + gamma_indices_.n_elem - 1).t() * cross_terms_).t();
    return all(d > (1.0 / n_));
}


double el_sem::get_objective()
{
    return -sum(log(d_));
}


void el_sem::set_gradient_hessian(vec &grad, mat &hessian, double b)
{
    mat cross_terms_divided = cross_terms_;
    cross_terms_divided.each_row() /= d_.t();
//    mat penalty_term_gamma = cross_terms_;
//    penalty_term_gamma.each_row() /= (b * (d_.t() - 1.0 / n_));

    mat residuals_divided_by_d = residuals_;
    residuals_divided_by_d.each_row() /= d_.t();
//    mat penalty_term_lambda = residuals_;
//    penalty_term_lambda.each_row() /= (b * (d_.t() - 1.0 / n_));

    grad.rows(0, v_ - 1) = -sum(residuals_divided_by_d, 1);// - penalty_term_lambda, 1);
    grad.rows(v_, v_ + gamma_indices_.n_elem - 1) = -sum(cross_terms_divided,1 ); // - penalty_term_gamma, 1);


    hessian.submat(0, 0, v_ - 1, v_ - 1) = residuals_divided_by_d * residuals_divided_by_d.t(); //+ penalty_term_lambda * penalty_term_lambda.t();
    hessian.submat(0, v_ , v_ - 1,  v_ + gamma_indices_.n_elem - 1) = residuals_divided_by_d * cross_terms_divided.t(); // + penalty_term_lambda * penalty_term_gamma.t() ;
    hessian.submat(v_, 0, v_ + gamma_indices_.n_elem - 1, v_ - 1) =  hessian.submat(0, v_ , v_ - 1,  v_ + gamma_indices_.n_elem - 1).t();
    hessian.submat(v_,  v_, v_ + gamma_indices_.n_elem -1, v_ + gamma_indices_.n_elem - 1) = cross_terms_divided * cross_terms_divided.t();// + penalty_term_gamma * penalty_term_gamma.t();

}

vec el_sem::get_d()
{
    return d_;
}
