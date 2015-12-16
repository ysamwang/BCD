#include "elSEMTies.h"

el_sem_tie:: el_sem_tie(SEXP b_weights_r, SEXP y_r, SEXP counts_r, SEXP omega_r,SEXP b_r , SEXP dual_r, int meanEst) : el_sem::el_sem(b_weights_r, y_r, omega_r, b_r , dual_r, meanEst)
{


    //matrix which holds the observed data
    counts_ = as<arma::vec>(counts_r);
    n_ = accu(counts_); //number of observations
}

double el_sem_tie::update_dual(double tol, int max_iter)
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
        el_sem_tie::set_gradient_hessian(grad, hessian);

        // back tracking line search

        update = solve(hessian, grad);

        back_tracking_counter = 0;
        while(!el_sem_tie::backtracking(update) && back_tracking_counter < max_back_track) {
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

    //constraints is r x q; dual is r x 1 so d is q x 1
    d_ = (constraints_.t() * dual_) + 1.0;
    if(conv_crit_ < tol) {
        return -sum(log(n_ * d_ /  counts_));
    } else {
        return -99999999;
    }

}

int el_sem_tie::backtracking(vec update)
{
    vec d = (constraints_.t() * (dual_ - update)) + 1.0;
    return all( d > (counts_ / n_));
}

void el_sem_tie::set_gradient_hessian(vec &grad, mat &hessian)
{
    int i;
    hessian.zeros();
    d_ = (constraints_.t() * dual_) + 1.0;
    grad = - sum( constraints_ * diagmat( counts_/  d_), 1);

    for(i = 0; i < counts_.n_elem ; i ++) {
        hessian += counts_(i) * constraints_.col(i) * constraints_.col(i).t() / pow(d_(i), 2);
    }
}

//return scaled d
vec el_sem_tie::get_d()
{
    return d_ * n_ / counts_;
}
