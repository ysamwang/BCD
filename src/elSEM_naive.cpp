#include "elSEM_naive.h"

el_sem_naive::el_sem_naive(SEXP weights_r, SEXP y_r, SEXP omega_r, SEXP b_r , SEXP dual_r, int meanEst)
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

    mat omega_weights(v_, v_, fill::zeros);
    mat b_weights(v_, v_, fill::zeros);
    vec means = vec(v_, fill::zeros); //non-zero means



    if(meanEst){
        means = as<arma::vec>(weights_r).head(v_);
        b_weights.elem(b_spots) = as<arma::vec>(weights_r).subvec(v_, v_ + b_spots.n_elem - 1);
        omega_weights.elem(omega_spots) = weights.subvec(v_ + b_spots.n_elem, v_ + b_spots.n_elem + omega_spots.n_elem - 1);
        omega_weights = symmatl(omega_weights);
    } else {
        b_weights.elem(b_spots) = weights.subvec(0, b_spots.n_elem - 1);
        omega_weights.elem(omega_spots) = weights.subvec(b_spots.n_elem, b_spots.n_elem + omega_spots.n_elem - 1);
        omega_weights = symmatl(omega_weights);
    }


    sigma_ = solve( (eye(v_, v_) - b_weights), solve( eye(v_, v_) - b_weights, omega_weights).t());
    dual_ = as<arma::vec>(dual_r);

    //
    constraints_ = mat(v_ + (v_ *(v_ + 1) )/ 2, n_);

    constraints_.rows(0, v_ - 1) = (eye(v_, v_) - b_weights) * y;
    if(meanEst) {
        constraints_.rows(0, v_ - 1).each_col() -= means;
    }

    int i,j,k;
    k = v_;
    for(i = 0; i < v_; i++) {
        for(j = 0; j <= i; j++ ) {
            constraints_.row(k) =  (y.row(i) % y.row(j)) - sigma_(i,j);
            k++;
        }
    }
    d_ = (constraints_.t() * dual_) + 1.0;
}

double el_sem_naive::update_dual(double tol, int max_iter)
{

    //pre-allocate memory for gradient, hessian and update
    vec grad(v_ + (v_ *(v_ + 1) )/ 2 );
    mat hessian(v_ + (v_ *(v_ + 1) )/ 2 ,v_ + (v_ *(v_ + 1) )/ 2 );
    vec update(v_ + (v_ *(v_ + 1) )/ 2 );
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

        //back track to ensure everything is greater than 1/n
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
        return -999999999;
    }

}

// check that update preserves all(d_ > 1/n)
int el_sem_naive::backtracking(vec update)
{
    vec d = (constraints_.t() * (dual_ - update)) + 1.0;
    return all( d > (1.0 / n_));
}


//update gradient and hessian
void el_sem_naive::set_gradient_hessian(vec &grad, mat &hessian)
{
    int i;
    hessian.zeros();
    d_ = (constraints_.t() * dual_) + 1.0;
    grad = -sum( constraints_ * diagmat( 1.0 / d_), 1);

    for(i = 0; i < n_ ; i ++) {
        hessian += constraints_.col(i) * constraints_.col(i).t() / pow(d_(i),2);
    }
}


//return scaled d
vec el_sem_naive::get_d()
{
    return d_ * n_;
}

//return convergence criteria
double el_sem_naive::get_conv_crit()
{
    return conv_crit_;
}

//return iterations used
int el_sem_naive::get_iter()
{
    return counter_;
}
