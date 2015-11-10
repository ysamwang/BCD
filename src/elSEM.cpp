#include "elSEM.h"

el_sem::el_sem(SEXP b_weights_r, SEXP y_r, SEXP omega_r, SEXP b_r , SEXP dual_r, int meanEst)
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
    b_weights_ = mat(v_, v_, fill::zeros); // matrix with edge weights


    if(meanEst){
        means = as<arma::vec>(b_weights_r).head(v_);
        b_weights_.elem(b_spots) = as<arma::vec>(b_weights_r).tail(as<arma::vec>(b_weights_r).n_elem - v_);
    } else {
    b_weights_.elem(b_spots) = as<arma::vec>(b_weights_r);
    }


    dual_ = as<arma::vec>(dual_r); //dual variables

    gamma_indices_ = arma::find(trimatu(as<arma::mat>(omega_r) == 0 )); //structural 0's in Omega

    constraints_ = mat(v_ + gamma_indices_.n_elem, n_);  //constraints containing mean and covariance restrictions
    constraints_.rows(0, v_ - 1) = (eye(v_, v_) - b_weights_) * y; //mean restrictions

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

double el_sem::update_dual(double tol, int max_iter)
{

    //pre-allocate memory for gradient, hessian and update
    vec grad(v_ + gamma_indices_.n_elem );
    mat hessian(v_ + gamma_indices_.n_elem , v_ + gamma_indices_.n_elem );
    vec update(v_ + gamma_indices_.n_elem );
    grad.zeros();
    hessian.zeros();
    update.zeros();


    // backtracking parameters
    double backtracking_scaling_const = .4;
    int back_tracking_counter;
    int max_back_track = 20; // steps required to scale update by 1e-8

    while(conv_crit_ > tol && counter_ < max_iter) {


        // build in backtracking if necessary
        set_gradient_hessian(grad, hessian);

        // back tracking line search

        update = solve(hessian, grad);

        back_tracking_counter = 0;
        while(!backtracking(update) && back_tracking_counter < max_back_track) {
            update *= backtracking_scaling_const;
            back_tracking_counter++;
        }
        //if back tracking does not terminate because of max iterations update
        //else terminate
        if(back_tracking_counter < max_back_track) {
            dual_ -= update;
        } else {
            return -99999;
        }


        conv_crit_ = norm(grad ,2);
        counter_++;
    }

    d_ = (constraints_.t() * dual_) + 1.0;
    if(conv_crit_ < tol) {
        return -sum(log(n_ * d_));
    } else {
        return -99999;
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
    grad = - sum( constraints_ * diagmat( 1.0 / d_), 1);

    for(i = 0; i < n_ ; i ++) {
        hessian += constraints_.col(i) * constraints_.col(i).t() / pow(d_(i),2);
    }
}

//return scaled d
vec el_sem::get_d()
{
    return d_ * n_;
}

//return convergence criteria
double el_sem::get_conv_crit()
{
    return conv_crit_;
}

//return iterations used
int el_sem::get_iter()
{
    return counter_;
}
