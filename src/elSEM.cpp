#include "elSEM.h"

el_sem::el_sem(SEXP y_r, SEXP omega_r, SEXP b_weights_r, SEXP dual_r, int v_r)
{
    y_ = as<arma::mat>(y_r);
    n_ = y_.n_cols;
    omega_ = as<arma::mat>(omega_r);
    b_weights_ = as<arma::mat>(b_weights_r);
    dual_ = as<arma::vec>(dual_r);
    v_ = v_r;

    gamma_indices_ = arma::find(trimatu(omega_ == 0 ));
    residuals_ = (eye(v_, v_) - b_weights_) * y_ ;
    cross_terms_.zeros(gamma_indices_.n_elem, n_);

    int i,j,k;
    for(k = 0; k < gamma_indices_.n_elem; k++) {
        j = (int) gamma_indices_(k) / v_;
        i = (int) gamma_indices_(k) % v_;
        cross_terms_.row(k) = residuals_.row(i) % residuals_.row(j);
    }
    constraints_ = join_cols(residuals_, cross_terms_);
    d_ = (constraints_.t() * dual_) + 1.0;
}

double el_sem::update_dual(double tol, int max_iter)
{
    double conv_crit = 1.0;
    int counter = 0;
    double backtracking_scaling_const = .7;

    //pre-allocate memory for gradient, hessian and update
    vec grad(v_ + gamma_indices_.n_elem );
    mat hessian(v_ + gamma_indices_.n_elem , v_ + gamma_indices_.n_elem );
    vec update(v_ + gamma_indices_.n_elem );
    grad.zeros();
    hessian.zeros();
    update.zeros();

    while(conv_crit > tol && counter < max_iter) {


        // build in backtracking if necessary
        set_gradient_hessian(grad, hessian);

//        Rcout << "Hessian: " <<endl << hessian <<endl;
//        Rcout << "Grad: " <<endl <<grad <<endl;
//        Rcout << cond(hessian)<< std::endl;


        // back tracking line search

        update = solve(hessian, grad);
        while(!backtracking(update)){
            update *= backtracking_scaling_const;
        }
        dual_ = dual_ - update;

        conv_crit = norm(grad ,2);
        counter++;
    }

    d_ = (constraints_.t() * dual_) + 1.0;
    if(conv_crit < tol){
        return -sum(log(n_ * d_));
    } else {
        return -9999;
    }

}

int el_sem::backtracking(vec update)
{
  vec d = (constraints_.t() * (dual_ + update)) + 1.0;
  return all( d > (1.0 / n_));
}

void el_sem::set_gradient_hessian(vec &grad, mat &hessian)
{
    int i;
    grad.zeros();
    hessian.zeros();
    for(i = 0; i < n_ ; i ++)
    {
        d_.row(i) = 1.0 + (dual_.t() * constraints_.col(i));
        grad -=  constraints_.col(i) / d_(i);
        hessian += constraints_.col(i) * constraints_.col(i).t() / pow(d_(i),2);
    }
}

vec el_sem::get_d()
{
    return d_ * n_;
}
