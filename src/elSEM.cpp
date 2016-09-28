#include "elSEM.h"


el_sem::el_sem(SEXP b_weights_r, SEXP y_r, SEXP b_r, bool mean_est_r, SEXP covar_restrict)
{
    y_ = as<arma::mat>(y_r); //matrix which holds the observed data
    n_ = y_.n_cols; //number of observations
    v_ = y_.n_rows; //number of variables
    counter_ = 0; //number of iterations
    conv_crit_ = 1.0; //convergence critiera for inner maximization (norm of gradient)
    mean_est_ = mean_est_r;

//    Rcpp::Rcout << "Point 1" << std::endl;
    //find appropriate spots to put in b_weights_r
    b_spots_ = find(as<arma::mat>(b_r)); //non-structural zeros in B
    vec means = vec(v_, fill::zeros); //non-zero means
    mat b_weights = mat(v_, v_, fill::zeros); // matrix with edge weights
    restrict_ = as<arma::mat>(covar_restrict); // matrix holding restrictions
//    Rcpp::Rcout << "Point 2" << std::endl;

    if(mean_est_) {
        means = as<arma::vec>(b_weights_r).head(v_);
        if(b_spots_.n_elem > 0){
            b_weights.elem(b_spots_) = as<arma::vec>(b_weights_r).tail(as<arma::vec>(b_weights_r).n_elem - v_);
        }
    } else {
        if(b_spots_.n_elem > 0){
            b_weights.elem(b_spots_) = as<arma::vec>(b_weights_r);
        }
    }
//        Rcpp::Rcout << "Point 3" << std::endl;
    // Filling in the constraint matrix
    constraints_ = mat(v_ + restrict_.n_rows, n_, fill::ones);
    // Mean restrictions
    constraints_.rows(0, v_ - 1) = (eye(v_, v_) - b_weights) * y_;
    // if there is a non-zero mean vector subtract off
    if(mean_est_) {
        constraints_.rows(0, v_ - 1).each_col() -= means;
    }
//    Rcpp::Rcout << "Point 4" << std::endl;

    // additional "uncorrelated" constraints
    // This could include covariance and arbitrary higher order constraints
    int k, z;
    for(k = 0; k < restrict_.n_rows; k++) {
        for(z = 0; z < restrict_.n_cols; z++) {
            // for lower order moments un-used elements of restrict_ will be -1
            // so only multiply in if not -1
            if(restrict_(k,z) != -1) {
                constraints_.row(k + v_) %= constraints_.row( restrict_(k, z) );
            }
        }
    }
//        Rcpp::Rcout << "Point 5" << std::endl;


    dual_ = vec(constraints_.n_rows, fill::zeros);

    // calculate denominator given dual variables
    // d_ = (constraints_.t() * dual_) + 1.0;
    d_ = vec(n_, fill::ones);
}

double el_sem::update_dual(double tol, int max_iter)
{
//            Rcpp::Rcout << "Point 6" << std::endl;

    vec grad(constraints_.n_rows, fill::zeros);
    vec update(constraints_.n_rows, fill::zeros);
    mat hessian(constraints_.n_rows, constraints_.n_rows, fill::zeros);

    // backtracking parameters
    int back_tracking_counter;
    int max_back_track = 20; // steps required to scale update by 1e-8


    while(conv_crit_ > tol) {

        // build in backtracking if necessary
        set_gradient_hessian(grad, hessian);

        update = solve(hessian, grad);

        back_tracking_counter = 0;
        while(!backtracking(update) && back_tracking_counter < max_back_track) {
            update *= BACKTRACKING_SCALING_CONST;
            if(++back_tracking_counter > max_back_track) {
                return INFEASIBLE_RETURN;
            }
        }
        // if back tracking does not return because of max back tracking iterations
        // update the dual
        dual_ -= update;

        // update the convergence criteria
        conv_crit_ = norm(grad ,2);

        // if max iterations hit, give infeasible return
        // In practice, this shouldn't be used often because
        // we have a strictly convex problem, but condition is built in for error checking
        if(counter_++ > max_iter) {
            return INFEASIBLE_RETURN;
        }
    }

    d_ = (constraints_.t() * dual_) + 1.0;
//            Rcpp::Rcout << "Point 7" << std::endl;

    return -sum(log(n_ * d_));
}


// check whether update remains in feasible space
int el_sem::backtracking(vec update)
{
    vec d = (constraints_.t() * (dual_ - update)) + 1.0;
    return all( d > (1.0 / n_));
}

// set gradient and hessian (wrt to dual variables) to appropriate values
void el_sem::set_gradient_hessian(vec &grad, mat &hessian)
{

    d_ = (constraints_.t() * dual_) + 1.0;
    grad = -sum( constraints_ * diagmat( 1.0 / d_), 1);

    hessian.zeros();
    int i;
    for(i = 0; i < n_ ; i ++) {
        hessian += constraints_.col(i) * constraints_.col(i).t() / pow(d_(i), 2);
    }
}

//return scaled d
vec el_sem::get_d()
{
    return d_ * n_;
}

vec el_sem::getGradient()
{

    // jacobian of the constraints wrt to the parameters
    mat dg_dtheta, dg_dtheta_n;

    // check if need to take derivative wrt to mean params as well
    if(mean_est_) {
        dg_dtheta = mat(constraints_.n_rows, v_ + b_spots_.n_elem);
        dg_dtheta_n = mat(constraints_.n_rows, v_ + b_spots_.n_elem);


        dg_dtheta.zeros();
        // loop through each observation
        int n, g, r, i, j, v, b;
        for(n = 0; n < n_; n++) {
            dg_dtheta_n.zeros();

            int v;
            // wrt to the mean parameters (just negative identity)
            for(v = 0; v < v_; v++) {
                dg_dtheta_n(v, v) = -1.0;
            }

            //loop through each free element of B
            for(b = 0; b < b_spots_.n_elem; b++) {
                j = b_spots_(b) / v_;
                i = b_spots_(b) % v_;

                // constraint i is the mean constraint for node i
                dg_dtheta_n(i, v_ + b) += -y_(j, n);

                // loop through uncorrelated constraints
                for(g = 0; g < restrict_.n_rows; g++) {
                    //loop through elements involved in constraints

                    for(r = 0; r < restrict_.n_cols; r++) {
                        // if one of the constraints involves epsilon_i
                        if(restrict_(g, r) == i) {
                            // take derivative wrt to beta_ij according to chain rule
                            // have effect of y * overall constraint where we've taken off the non-product terms
                            dg_dtheta_n(v_ + g, v_ + b) += -y_(j, n) * constraints_(v_ + g, n) / constraints_(i, n);
                        }
                    }
                }
            }

            // sum up
            dg_dtheta += dg_dtheta_n / d_(n);
        }
        // end meanEst
    } else { // no mean params

        // set up dg_dtheta holder
        dg_dtheta = mat(constraints_.n_rows, b_spots_.n_elem);
        dg_dtheta.zeros();

        dg_dtheta_n = mat(constraints_.n_rows, b_spots_.n_elem);

        int n, g, r, i, j, b;
        // loop through each observation
        for(n = 0; n < n_; n++) {
            dg_dtheta_n.zeros(); //reset to 0

            //loop through each non-structural zero of B
            for(b = 0; b < b_spots_.n_elem; b++) {
                // row (i) and column (j) of particular element of b
                j = b_spots_(b) / v_;
                i = b_spots_(b) % v_;

                // constraint i is the mean constraint for node i
                dg_dtheta_n(i, b) += -y_(j, n);

                // loop through uncorrelated constraints
                for(g = 0; g < restrict_.n_rows; g++) {
                    //loop through elements involved in constraints

                    for(r = 0; r < restrict_.n_cols; r++) {
                        // if one of the constraints involves epsilon_i
                        if(restrict_(g, r) == i) {
                            // take derivative wrt to beta_ij according to chain rule
                            // have effect of y * overall constraint where we've taken off the non-product terms
                            dg_dtheta_n(v_ + g, b) += -y_(j, n) * constraints_(v_ + g, n) / constraints_(i, n);
                        }
                    }
                }
            }

            // sum up by weight
            dg_dtheta += dg_dtheta_n / d_(n);
        }
    }

//            Rcpp::Rcout << "Point 9" << std::endl;

//    Rcpp::Rcout << dg_dtheta <<std::endl;
//    Rcpp::Rcout << dual_ <<std::endl;
    mat logEL_grad_wrt_b = -(dual_.t() * dg_dtheta) ;
//    Rcpp::Rcout << "Point 9.5" << std::endl;

    return logEL_grad_wrt_b.row(0).t();
}

