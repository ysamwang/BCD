#include "elSEM.h"

el_sem::el_sem() {}

el_sem::el_sem(SEXP b_weights_r, SEXP y_r, SEXP b_r, SEXP moment_2_restrictions_r, SEXP moment_3_restrictions_r,  SEXP moment_4_restrictions_r, int meanEst)
{
    mat y = as<arma::mat>(y_r); //matrix which holds the observed data
    n_ = y.n_cols; //number of observations
    v_ = y.n_rows; //number of variables
    counter_ = 0; //number of iterations
    conv_crit_ = 1.0; //convergence critiera (norm of gradient)

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

    dual_ = vec(constraints_.n_rows, fill::zeros);


    // Filling in the constraint matrix

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

    d_ = (constraints_.t() * dual_) + 1.0;
}

double el_sem::update_dual(double tol, int max_iter)
{
    vec grad(constraints_.n_rows, fill::zeros);
    vec update(constraints_.n_rows, fill::zeros);
    mat hessian(constraints_.n_rows, constraints_.n_rows, fill::zeros);

    // backtracking parameters
    double backtracking_scaling_const = .4;
    int back_tracking_counter;
    int max_back_track = 20; // steps required to scale update by 1e-8

    while(conv_crit_ > tol && counter_ < max_iter) {


        // build in backtracking if necessary
        set_gradient_hessian(grad, hessian);

        // back tracking line search
        if (cond(hessian) > 1e4) {
//            Rcout << "d: " << endl << d_ << endl;
        }

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
        hessian += constraints_.col(i) * constraints_.col(i).t() / pow(d_(i), 2);
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


//return iterations used
vec el_sem::get_dual()
{
    return dual_;
}

vec el_sem::getGradient()
{
  vec logEL_gradient(dual_.n_rows);
  logEL_gradient.zeros();
  int n;
  
  for(n = 0; n < n_; n++){
    logEL_gradient += dual_.t() *  /d_(n);
  }
  
  
  //need to finish
  
}