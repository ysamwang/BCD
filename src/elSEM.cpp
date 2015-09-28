#include "elSEM.h"

el_sem::el_sem(SEXP y_r, SEXP b_r, SEXP omega_r, SEXP b_weights_r, SEXP omega_weights_r, SEXP d_r, SEXP lambda_r, SEXP gamma_r, int v_r)
{
    y_ = as<arma::mat>(y_r);
    n_ = y_.n_cols;
    b_ = as<arma::mat>(b_r);
    omega_ = as<arma::mat>(omega_r);
    b_weights_ = as<arma::mat>(b_weights_r);
    omega_weights_ = as<arma::mat>(omega_weights_r);
    d_ = as<arma::vec>(d_r);
    lambda_ = as<arma::vec>(lambda_r);
    gamma_ = as<arma::vec>(gamma_r);
    v_ = v_r;

    gamma_indices_ = arma::find(omega_ == 0);
    residuals_ = (eye(v_, v_) - b_weights_) * y_ ;
    d_.zeros(n_);

    cross_terms_.zeros(accu(omega_ == 0), n_);

    int i,j,k;
    for(k = 0; k < gamma_indices_.n_elem; k++) {
        j = (int) gamma_indices_(k) / v_;
        i = (int) gamma_indices_(k) % v_;
        cross_terms_.row(k) = residuals_.row(i) % residuals_.row(j);
    }
    Rcout <<" Initialization Done" <<std::endl;

}

void el_sem::update_lambda_gamma(double tol, int max_iter)
{
    double conv_crit = 1.0;
    int counter = 0;
    double old_obj, new_obj;
    new_obj = get_objective();

    //pre-allocate memory for gradient, hessian and update
    vec grad(v_ + gamma_indices_.n_elem);
    mat hessian(v_ + gamma_indices_.n_elem, v_ + gamma_indices_.n_elem);
    vec update(v_ + gamma_indices_.n_elem);
    set_d_star();
    while(conv_crit > tol && counter < max_iter) {
        old_obj = new_obj;


        // build in backtracking if necessary
        set_gradient_hessian(grad, hessian);

        Rcout <<"Grad Set" <<std::endl;
        update =  join_cols(lambda_, gamma_) - solve(hessian, grad);


        lambda_ = update.rows(0, v_ - 1);
        gamma_ = update.rows(v_, v_ + gamma_indices_.n_elem - 1);
        Rcout << "Updated Lambda and Gamma" <<std::endl;
        set_d_star();
        Rcout << "Updated d" <<std::endl;
        new_obj = get_objective();
        conv_crit = fabs((new_obj - old_obj)/old_obj);
        counter++;
        Rcout <<"Loop Done: " << new_obj <<endl;
    }
}

void el_sem::set_d_star()
{
    d_.ones();
    // term from


    Rcout << lambda_ <<std::endl;
    Rcout << residuals_ <<std::endl;
    d_ += (lambda_.t() * residuals_).t();
    Rcout << "Update Lambda" <<std::endl;
    d_ += (gamma_.t() * cross_terms_).t();
    Rcout << "Update Gamma" <<std::endl;
 }



double el_sem::get_objective()
{
    return -sum(log(d_));
}


void el_sem::set_gradient_hessian(vec grad, mat hessian)
{
    mat cross_terms_divided = cross_terms_;
    cross_terms_divided.each_row() /= d_.t();
    mat residuals_divided_by_d = residuals_;
    residuals_divided_by_d.each_row() /= d_.t();

    Rcout <<"Updating First Step" <<std::endl;
    grad.rows(0, v_ - 1) = sum(residuals_divided_by_d, 1);
    grad.rows(v_, v_ + gamma_indices_.n_elem - 1) = sum(cross_terms_divided, 1);
    Rcout <<"Updating Grad" <<std::endl;

    hessian.submat(0, 0, v_ - 1, v_ - 1) = residuals_divided_by_d * residuals_divided_by_d.t();
    hessian.submat(0, v_ , v_ - 1,  v_ + gamma_indices_.n_elem - 1) = residuals_divided_by_d * cross_terms_divided.t();
    hessian.submat(v_, 0, v_ + gamma_indices_.n_elem - 1, v_ - 1) =  hessian.submat(0, v_ , v_ - 1,  v_ + gamma_indices_.n_elem - 1).t();
    hessian.submat(v_,  v_, v_ + gamma_indices_.n_elem -1, v_ + gamma_indices_.n_elem - 1) = cross_terms_divided * cross_terms_divided.t();
    Rcout <<"Updating Hessian" <<std::endl;
}
