#include "elSEM_naive.h"

el_sem_naive::el_sem_naive(SEXP weights_r, SEXP y_r, SEXP omega_r, SEXP b_r , SEXP dual_r)
{
    mat y = as<arma::mat>(y_r);
    n_ = y.n_cols;
    v_ = y.n_rows;
    counter_ = 0;
    conv_crit_ = 1.0;

    //find appropriate spots to put in b_weights_r and omega_weights so that we can build sigma
    uvec b_spots = find(as<arma::mat>(b_r));
    uvec omega_spots = arma::find(trimatl(as<arma::mat>(omega_r)));
    vec weights = as<arma::vec>(weights_r);

    mat b_weights(v_, v_, fill::zeros);
    b_weights.elem(b_spots) = weights.rows(0, b_spots.n_elem - 1);


    mat omega_weights(v_, v_, fill::zeros);
    omega_weights.elem(omega_spots) = weights.rows(b_spots.n_elem, b_spots.n_elem + omega_spots.n_elem -1);
    omega_weights = symmatl(omega_weights);
    sigma_ = solve( (eye(v_, v_) - b_weights), solve( eye(v_, v_) - b_weights, omega_weights).t());
    dual_ = as<arma::vec>(dual_r);

    //
    constraints_ = mat(v_ + (v_ *(v_ + 1) )/ 2, n_);

    constraints_.rows(0, v_ - 1) = (eye(v_, v_) - b_weights) * y;
    int i,j,k;
    k = v_;
    for(i = 0; i < v_; i++){
        for(j = 0; j <= i; j++ ) {
            constraints_.row(k) =  (y.row(i) % y.row(j)) - sigma_(i,j);
            k++;
        }
    }
    d_ = (constraints_.t() * dual_) + 1.0;
}

double el_sem_naive::update_dual(double tol, int max_iter)
{
    double backtracking_scaling_const = .7;

    //pre-allocate memory for gradient, hessian and update
    vec grad(v_ + (v_ *(v_ + 1) )/ 2 );
    mat hessian(v_ + (v_ *(v_ + 1) )/ 2 ,v_ + (v_ *(v_ + 1) )/ 2 );
    vec update(v_ + (v_ *(v_ + 1) )/ 2 );
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
        return -999999999;
    }

}

int el_sem_naive::backtracking(vec update)
{
  vec d = (constraints_.t() * (dual_ - update)) + 1.0;
  return all( d > (1.0 / n_));
}

void el_sem_naive::set_gradient_hessian(vec &grad, mat &hessian)
{
    int i;
    hessian.zeros();
    d_ = (constraints_.t() * dual_) + 1.0;
    grad = -sum( constraints_ * diagmat( 1.0 / d_), 1);

    for(i = 0; i < n_ ; i ++)
    {
        hessian += constraints_.col(i) * constraints_.col(i).t() / pow(d_(i),2);
    }
}

vec el_sem_naive::get_d()
{
    return d_ * n_;
}

double el_sem_naive::get_conv_crit()
{
    return conv_crit_;
}

int el_sem_naive::get_iter()
{
    return counter_;
}
