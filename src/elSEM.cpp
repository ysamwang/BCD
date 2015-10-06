#include "elSEM.h"

el_sem::el_sem(SEXP b_weights_r, SEXP y_r, SEXP omega_r, SEXP b_r , SEXP dual_r)
{
    mat y = as<arma::mat>(y_r);
    n_ = y.n_cols;
    v_ = y.n_rows;
    counter_ = 0;
    conv_crit_ = 1.0;

    //find appropriate spots to put in b_weights_r
    uvec b_spots = find(as<arma::mat>(b_r));
    b_weights_ = mat(v_, v_, fill::zeros);
    b_weights_.elem(b_spots) = as<arma::vec>(b_weights_r);

    dual_ = as<arma::vec>(dual_r);

    gamma_indices_ = arma::find(trimatu(as<arma::mat>(omega_r) == 0 ));
    constraints_ = mat(v_ + gamma_indices_.n_elem, n_);
    constraints_.rows(0, v_ - 1) = (eye(v_, v_) - b_weights_) * y ;

    int i,j,k;
    for(k = 0; k < gamma_indices_.n_elem; k++) {
        j = (int) gamma_indices_(k) / v_;
        i = (int) gamma_indices_(k) % v_;
        constraints_.row(k + v_) = constraints_.row(i) % constraints_.row(j);
    }
    d_ = (constraints_.t() * dual_) + 1.0;
}

double el_sem::update_dual(double tol, int max_iter)
{
    double backtracking_scaling_const = .7;

    //pre-allocate memory for gradient, hessian and update
    vec grad(v_ + gamma_indices_.n_elem );
    mat hessian(v_ + gamma_indices_.n_elem , v_ + gamma_indices_.n_elem );
    vec update(v_ + gamma_indices_.n_elem );
    grad.zeros();
    hessian.zeros();
    update.zeros();

    while(conv_crit_ > tol && counter_ < max_iter) {


        // build in backtracking if necessary
        set_gradient_hessian(grad, hessian);

        // back tracking line search

        update = solve(hessian, grad);
        while(!backtracking(update)){
            update *= backtracking_scaling_const;
        }
        dual_ -= update;

        conv_crit_ = norm(grad ,2);
        counter_++;
    }

    d_ = (constraints_.t() * dual_) + 1.0;
    if(conv_crit_ < tol){
        return -sum(log(n_ * d_));
    } else {
        return -9999;
    }

}

int el_sem::backtracking(vec update)
{
  vec d = (constraints_.t() * (dual_ - update)) + 1.0;
  return all( d > (1.0 / n_));
}

void el_sem::set_gradient_hessian(vec &grad, mat &hessian)
{
    int i;
    hessian.zeros();
    d_ = (constraints_.t() * dual_) + 1.0;
    grad = - sum( constraints_ *diagmat( 1.0 / d_), 1);

    for(i = 0; i < n_ ; i ++)
    {
        hessian += constraints_.col(i) * constraints_.col(i).t() / pow(d_(i),2);
    }
}

vec el_sem::get_d()
{
    return d_ * n_;
}

double el_sem::get_conv_crit()
{
    return conv_crit_;
}

int el_sem::get_iter()
{
    return counter_;
}
