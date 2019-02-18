#include "elSEM_naive.h"

el_sem_naive::el_sem_naive(SEXP weights_r, SEXP y_r, SEXP b_r, SEXP omega_r, bool mean_est_r)
{
    vec means;
    y_ = as<arma::mat>(y_r);
    n_ = y_.n_cols;
    v_ = y_.n_rows;
    mat b_weights(v_, v_, fill::zeros);
    mat omega_weights(v_, v_, fill::zeros);
    mean_est_ = mean_est_r;


    //find appropriate spots to put in b_weights_r and omega_weights so that we can build sigma
    b_spots_ = find(as<arma::mat>(b_r));
    omega_spots_ = arma::find(trimatl(as<arma::mat>(omega_r)));
    vec weights = as<arma::vec>(weights_r);

    if(mean_est_) {
        means = as<arma::vec>(weights_r).head(v_);
        b_weights.elem(b_spots_) = as<arma::vec>(weights_r).subvec(v_, v_ + b_spots_.n_elem - 1);
        omega_weights.elem(omega_spots_) = weights.subvec(v_ + b_spots_.n_elem, v_ + b_spots_.n_elem + omega_spots_.n_elem - 1);
        omega_weights = symmatl(omega_weights);
    } else {
        b_weights.elem(b_spots_) = weights.subvec(0, b_spots_.n_elem - 1);
        omega_weights.elem(omega_spots_) = weights.subvec(b_spots_.n_elem, b_spots_.n_elem + omega_spots_.n_elem - 1);
        omega_weights = symmatl(omega_weights);
    }

    // Create matrix of constraints
    // set mean constraints
    constraints_ = mat(v_ + (v_ *(v_ + 1) )/ 2, n_);
    constraints_.rows(0, v_ - 1) = (eye(v_, v_) - b_weights) * y_;

    if(mean_est_) {
        constraints_.rows(0, v_ - 1).each_col() -= means;
    }

    int i,j,k;
    k = v_;
    for(i = 0; i < v_; i++) {
        for(j = 0; j <= i; j++ ) {
            constraints_.row(k) =  (constraints_.row(i) % constraints_.row(j)) - omega_weights(i, j);
            k++;
        }
    }
    //dual variables
    dual_ = vec(constraints_.n_rows, fill::zeros);
//    d_ = (constraints_.t() * dual_) + 1.0;
    d_ = vec(n_, fill::ones);

}

double el_sem_naive::update_dual(double tol, int max_iter)
{

    //pre-allocate memory for gradient, hessian and update
    vec grad(constraints_.n_rows );
    vec update(constraints_.n_rows);
    mat hessian(constraints_.n_rows , constraints_.n_rows);
    grad.zeros();
    hessian.zeros();
    update.zeros();

    // backtracking parameters
    int back_tracking_counter;
    int max_back_track = 20; // steps required to scale update by 1e-8


    // search parameters
    double conv_crit = tol + 1.0;
    int counter = 0;

    while(conv_crit > tol) {

        // build in backtracking if necessary
        set_gradient_hessian(grad, hessian);

        // back tracking line search
        update = solve(hessian, grad);

        //back track to ensure everything is greater than 1/n
        back_tracking_counter = 0;
        while(!backtracking(update) ) {
            update *= BACKTRACKING_SCALING_CONST;
            if(++back_tracking_counter > max_back_track) {
                return INFEASIBLE_RETURN;
            }
        }

        dual_ -= update;

        conv_crit = norm(grad ,2);
        if(counter++ > max_iter) {
            return INFEASIBLE_RETURN;
        }
    }

    d_ = (constraints_.t() * dual_) + 1.0;
    return -sum(log(n_ * d_));
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


vec el_sem_naive::getGradient()
{
    // jacobian of the constraints wrt to the parameters
    mat dg_dtheta, dg_dtheta_n;
    // check if need to take derivative wrt to mean params as well
    if(mean_est_) {
        dg_dtheta = mat(constraints_.n_rows, v_ + b_spots_.n_elem + omega_spots_.n_elem);
        dg_dtheta_n = mat(constraints_.n_rows, v_ + b_spots_.n_elem + omega_spots_.n_elem);

        dg_dtheta.zeros();
        // loop through each observation
        int n, v, b, i, j, z, k, w;
        for(n = 0; n < n_; n++) {
            dg_dtheta_n.zeros();
            // wrt to the mean parameters (just negative identity)
            for(v = 0; v < v_; v++) {
                dg_dtheta_n(v, v) = -1.0;
            }

            //loop through each free element of B
            for(b = 0; b < b_spots_.n_elem; b++) {
                j = b_spots_(b) / v_;
                i = b_spots_(b) % v_;

                // mean constraint
                dg_dtheta_n(i,v_ +  b) += -y_(j, n);

                // hit all covariance constraints
                for(z = 0; z < v_; z++) {
                    if (z <= i) {
                        k = i * (i + 1) / 2 - (i - z) - 1;
                    } else  {
                        k = z* (z + 1) / 2 - (z - i) - 1;
                    }

                    if (z != i) {
                        dg_dtheta_n(v_ + k, v_ + b) += -y_(j, n) * constraints_(z, n);
                    } else {
                        dg_dtheta_n(v_ + k, v_ + b) += -2 * y_(j, n) * constraints_(z, n);
                    }
                }
            }

            // loop through each free element of Omega
            for(w = 0; w < omega_spots_.n_elem; w++) {
                j = omega_spots_(w) / v_;
                i = omega_spots_(w) % v_;

                if (j <= i) {
                    k = i * (i + 1) / 2 - (i - j) - 1;
                } else  {
                    k = j * (j + 1) / 2 - (j - i) - 1;
                }

                dg_dtheta_n(v_ + k, v_ + b_spots_.n_elem + w) += -1.0;
            }

            dg_dtheta += dg_dtheta_n / d_(n);
        }
        // end meanEst
    } else { // no mean params

        dg_dtheta = mat(constraints_.n_rows, b_spots_.n_elem + omega_spots_.n_elem);
        dg_dtheta_n = mat(constraints_.n_rows, b_spots_.n_elem + omega_spots_.n_elem);

        dg_dtheta.zeros();


        // loop through each observation
        int n, v, z, k, w, b, i, j;
        for(n = 0; n < n_; n++) {
            dg_dtheta_n.zeros();

            //loop through each free element of B
            for(b = 0; b < b_spots_.n_elem; b++) {
                j = b_spots_(b) / v_;
                i = b_spots_(b) % v_;

                // derivative of mean constraint wrt b_{ij}
                dg_dtheta_n(i, b) = -y_(j, n);

                // hit all covariance constraints
                for(z = 0; z < v_; z++) {
                    if (z <= i) {
                        k = (i+1) * (i + 2) / 2 - (i - z) - 1;
                    } else  {
                        k = (z + 1) * (z + 2) / 2 - (z - i) - 1;
                    }

                    if (z != i) {
                        dg_dtheta_n(v_ + k, b) += -y_(j, n) * constraints_(z, n);
                    } else {
                        dg_dtheta_n(v_ + k, b) += -2 * y_(j, n) * constraints_(z, n);
                    }
                }
            }

            // loop through each free element of Omega
            for(w = 0; w < omega_spots_.n_elem; w++) {
                j = omega_spots_(w) / v_;
                i = omega_spots_(w) % v_;


                if (j <= i) {
                    k = (i+1) * (i + 2) / 2 - (i - j) - 1;
                } else  {
                    k = (j + 1) * (j + 2) / 2 - (j - i) - 1;
                }

//                Rcout << i << " " << j  << " " << k <<std::endl;

                dg_dtheta_n(v_ + k, b_spots_.n_elem + w) += -1.0;
            }

            dg_dtheta += dg_dtheta_n / d_(n);
//            Rcout << dg_dtheta_n <<std::endl;
        }
        // end meanEst
    }

    mat logEL_grad_wrt_b = -(dual_.t() * dg_dtheta) ;

    return logEL_grad_wrt_b.row(0).t();
}

