#include "linSEM.h"

using namespace Rcpp ;
using namespace arma ;

//Constructor
// Br is a matrix of 1's an 0's indicating the directed edge structure
// Omegar is a matrix of 1's and 0's indicatin the bidirected edge structure
// BInitr is a matrix of initial directed edge weights
// OmegaInitr is a matrix of initial bidirected edge weights
// Yr is an V x n matrix of observations
lin_sem::lin_sem(SEXP Br, SEXP Omegar, SEXP BInitr, SEXP OmegaInitr, SEXP Yr, double OmegaInitScale)
{
    // take SEXP data structures and convert to armadillo data structures
    B = as<arma::mat>(Br);
    Omega =  as<arma::mat>(Omegar);
    Y = as<arma::mat>(Yr);
    V = B.n_cols; // number of nodes

    //if BInit is null, use initialization routine
    // else use given initialization
    if(Rf_isNull(BInitr)) {
        initEst(OmegaInitScale);
    } else {
        // clone since these will be updated
        BInit = as<arma::mat>(clone(BInitr));
        OmegaInit =  as<arma::mat>(clone(OmegaInitr));
    }
    //initialize singleUpdates to all 0's
    singleUpdates.zeros(V);
}

void lin_sem::initEst(double OmegaInitScale)
{
    int i;

    //initialize to all 0's with
    BInit.zeros(V, V);
    OmegaInit.zeros(V,V);

    //initialize B and Omega
    uvec i_vec(1); // uvec to get i'th row/column
    mat resid; //matrix of residuals from regressions
    resid.copy_size(Y); //same size as Y

    for(i = 0; i < V; i++) {
        // Set B through OLS
        uvec pa_i = find(B.row(i));
        if( pa_i.n_rows > 0) {
            vec coef = solve(Y.rows(pa_i).t(), Y.row(i).t());
            i_vec[0] = i;
            BInit(i_vec, pa_i ) = coef.t();

            //Estimate residual from OLS and use as initial estimates
            resid.row(i) = (Y.row(i) - coef.t() * Y.rows(pa_i));
        } else {
            //if there are no parents
            resid.row(i) = Y.row(i);
        }
    }

    // punch 0's into the sample covariance and get eigen values
    OmegaInit = Omega % cov(resid.t());
    vec eigval = eig_sym( OmegaInit );


    //Check if OmegaInit is PD
    if(!all(eigval > 0)) {
        //if not, scale down rows so that it is diagonally dominant to ensure PD
        rowvec row_i;
        double row_sum;
        int j;

        for(i = 0; i < V; i++) {
            //get sum of row i without element i
            row_i = OmegaInit.row(i);
            row_i.shed_col(i);
            row_sum = sum(abs(row_i));

            //check for diagonal dominance in row i
            if(row_sum > OmegaInit(i,i)) {
                //scale down each element so that the sum of the off-diagonals is equal to (omega_ii * OmegaInitScale)
                for(j = 0; j < V; j++) {
                    if(j != i) {
                        OmegaInit(i, j) = OmegaInit(i, j) * OmegaInit(i,i) * OmegaInitScale / row_sum;
                        OmegaInit(j, i) = OmegaInit(i, j);
                    }
                }
            }
        }
    }
}


// Main function for lin_sem object
// Updates relevant parameters for node i
// takes node i and max condition number
// Main function for lin_sem object
// Updates relevant parameters for node i
// takes node i and max condition number
int lin_sem::updateNode(int i, double maxKappa)
{
    // get parents and siblings
    uvec pa_i = find(B.row(i));
    uvec sib_i = find(Omega.row(i));
    arma::uvec i_in_sib_i = find(sib_i == i);
    int i_index = i_in_sib_i[0];
    sib_i.shed_row(i_index); //remove self from siblings

    int s = sib_i.n_elem; //number of siblings and parents
    int p = pa_i.n_elem;

    //uvec which hold's i
    uvec i_vec(1);
    i_vec[0] =  i;

    arma::mat X;
    arma::mat Z;
    arma::vec aHat, solved;
    arma::rowvec resid;
    arma::uvec sib2;
    double w_ii_two_norm;
    int k;

    // Check conditioning number
    double condNumber = cond(omegaSubset(i));
    if(condNumber > maxKappa) {
        Rcout <<"Condition Number too large: " << condNumber <<endl;
        return 0;
    }




    //calculate c_0 and c_pa
    arma::vec c_pa(pa_i.n_elem, fill::zeros);

    arma::mat subMat_cPa;
    for(k = 0; k < pa_i.n_elem; k ++) {
        subMat_cPa = eyeBSubset(i);
        subMat_cPa.shed_col( pa_i(k));
        c_pa(k) = det(subMat_cPa) * pow( -1.0, i + pa_i(k) + 1);
    }

    BInit.row(i) = vec(V, fill::zeros).t();
    arma::mat subMat_cNaught = eyeBSubset(i);
    subMat_cNaught.shed_col(i);

    double cNaught = det(subMat_cNaught);

    // Form X_i matrix
    if(s > 0 ) {

        //get psuedo variables
        arma::uvec sib2 = sib_i;
        sib2.elem(find(sib2 > i)) = sib2.elem(find(sib2 > i)) - 1;

        arma::mat Z = solve(omegaSubset(i), eyeBSubset(i) * Y);
        if( p > 0 ) {
            // parents and siblings
            // combine into a single matrix
            X = join_vert(Z.rows(sib2), Y.rows(pa_i));
        } else {
            //siblings no parents
            X = Z.rows(sib2);
        }
    } else {
        // no siblings but parents
        if(p > 0) {
            X = Y.rows(pa_i);
        } else {
            resid = Y.row(i);
        }
    }





    //Solve for minimizer of all parameters
    if(s + p > 0) {
        aHat = solve(X.t(), Y.row(i).t());
        resid = Y.row(i) - aHat.t() * X;


        if(norm(c_pa) == 0) {
            solved = aHat;

            // if there are no directed cycles (which we assume is true if norm(a_pa) == 0)
            // and no siblings, then you only need one update
            if(s == 0) {
//                singleUpdates(i) = 1;
            }
        } else {

            vec augmented_cPa = join_cols(vec(s, fill::zeros), c_pa);

            vec back_half = pow(norm(resid), 2) / (cNaught + dot(augmented_cPa, aHat)) * solve(X * X.t(), augmented_cPa);
            solved = aHat + back_half;
        }

        w_ii_two_norm = pow(norm(Y.row(i) - solved.t() * X),2);

    } else {
        w_ii_two_norm = pow(norm(Y.row(i)),2);
    }




    //Update the parameters
    if(p > 0) {
        BInit(i_vec, pa_i) = solved.rows(s, s + p - 1).t();
    }

    if(s > 0) {
        OmegaInit(i_vec, sib_i) = solved.rows(0, s - 1).t();
        OmegaInit(sib_i, i_vec) = solved.rows(0, s - 1);
    }

    // Update variance term
    arma::vec omegaVec = OmegaInit.col(i);
    omegaVec.shed_row(i);


        OmegaInit(i, i) = w_ii_two_norm / Y.n_cols +
                          dot(omegaVec.t() , solve(omegaSubset(i),  omegaVec));




    // succesfully updated. If condition number of Omega[-i,-i] is too large, then it returns 0 and does not update
    return 1;
}



mat lin_sem::omegaSubset(int i)
{
    arma::mat ret;
    ret = OmegaInit + 0.0;
    ret.shed_row(i);
    ret.shed_col(i);
    return ret;
}


mat lin_sem::eyeBSubset(int i)
{
    arma::mat ret = eye(V, V) - BInit;
    ret.shed_row(i);
    return ret;
}


mat lin_sem::getBInit()
{
    return BInit;
}

mat lin_sem::getOmegaInit()
{
    return OmegaInit;
}

mat lin_sem::getSigma()
{
    return solve(eye(V, V) - BInit, solve(eye(V, V) - BInit, OmegaInit).t());
}

int lin_sem::getV()
{
    return V;
}

int lin_sem::singleUpdateOnly(int i)
{
    return (int) singleUpdates[i];
}

mat lin_sem::getB()
{
    return B;
}

mat lin_sem::getOmega()
{
    return Omega;
}
