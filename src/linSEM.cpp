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
    if(Rf_isNull(BInitr))
    {
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
            //if there are no parents, just get deviations from mean
            resid.row(i) = Y.row(i) - mean(Y.row(i));
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

      for(i = 0; i < V; i++){
        //get sum of row i without element i
        row_i = OmegaInit.row(i);
        row_i.shed_col(i);
        row_sum = sum(abs(row_i));

        //check for diagonal dominance in row i
        if(row_sum > OmegaInit(i,i))
        {
          //scale down each element so that the sum of the off-diagonals is equal to (omega_ii * OmegaInitScale)
          for(j = 0; j < V; j++)
          {
            if(j != i){
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
    arma::vec solved;
    arma::rowvec resid;
    arma::uvec sib2;
    int k;
    double det1, det2, aNaught;

    double condNumber = cond(omegaSubset(i));
    if(condNumber > maxKappa){
      Rcout <<"Condition Number too large: " << condNumber <<endl;
      return 0;
    }


    //calculate a_0 and a_pa
    arma::vec a_pa(pa_i.n_elem, fill::zeros);

    for(k = 0; k < pa_i.n_elem; k ++) {
        //potential speed up by identifying strongly connected components
        BInit(i, pa_i(k)) = 2.0;
        det2 = det(eye(V,V) - BInit);
        BInit(i, pa_i(k)) = 1.0;
        det1 = det(eye(V,V) - BInit);
        a_pa(k) = det2 - det1;
    }

    BInit.row(i) = vec(V, fill::zeros).t();
    arma::mat subMat_aNaught = eye(V, V) - BInit;
    subMat_aNaught.shed_row(i);
    subMat_aNaught.shed_col(i);
    aNaught = det(subMat_aNaught);

    //Acyclic case or case where no parents are involved in cycles
    if(norm(a_pa, 2) == 0) {

       // There are some siblings
      if(s > 0 ){

        //get psuedo variables
        arma::uvec sib2 = sib_i;
        sib2.elem(find(sib2 > i)) = sib2.elem(find(sib2 > i)) - 1;

        arma::mat Z = solve(omegaSubset(i), eyeBSubset(i)) * Y;
        if( p > 0 ) {
        // combine into a single matrix
          X = join_vert(Z.rows(sib2), Y.rows(pa_i));
          solved = solve(X.t(), Y.row(i).t());

          //update Binit since there are parents
          BInit(i_vec, pa_i) = solved.rows(s, s + p - 1).t();
        } else {

          X = Z.rows(sib2);
          solved = solve(X.t(), Y.row(i).t());
        }

        //update OmegaInit since there are siblings
        OmegaInit(i_vec, sib_i) = solved.rows(0, s - 1).t();
        OmegaInit(sib_i, i_vec) = solved.rows(0, s - 1);

        resid = Y.row(i) - solved.t() * X;


      } else { //no siblings

        //since there are no siblings and a_pa = 0, then this only needs to be updated once
        singleUpdates(i) = 1;

        if( p > 0 ) {
          X = Y.rows(pa_i);
          solved = solve(X.t(), Y.row(i).t());
          BInit(i_vec, pa_i) = solved.rows(0, p - 1).t();
          resid = Y.row(i) - solved.t() * X;
        } else {
        // Assume mean 0, but not necessarily centered data
          resid = Y.row(i);
        }
      }

      arma::vec omegaVec = OmegaInit.col(i);
      omegaVec.shed_row(i);



      OmegaInit(i, i) = 1.0 / Y.n_cols * pow (norm(resid), 2) +
        dot(omegaVec.t() , solve(omegaSubset(i),  omegaVec));

    } else {
        // Full BCD
        //get Qp
        arma::mat Qp;

        //normalize a_pa and aNaught
        aNaught = aNaught / norm(a_pa, 2);
        a_pa = a_pa / norm(a_pa, 2);


        if(a_pa.n_elem > 1) {
            Qp =  join_vert(null(a_pa.t()).t(), a_pa.t());
        } else {
          Qp << a_pa(0) <<endr;
        }


        if(s > 0){
          //get psuedo variables
          sib2 = sib_i;
          sib2.elem(find(sib2 > i)) = sib2.elem(find(sib2 > i)) - 1;
          Z = solve(omegaSubset(i), eyeBSubset(i)) * Y;
          X = join_vert(Z.rows(sib2), Qp * Y.rows(pa_i));

        } else {
          X =  Qp * Y.rows(pa_i);
        }


        // QR decomposition of X^t
        arma::mat Q, R;

        qr(Q, R, X.t());

                // Compute y_0^2 for gamma^\star
        double yNaught2 = 0;
        for(k = s + p; k < Y.n_cols; k++) {
            yNaught2 += pow(dot(Y.row(i), Q.col(k)), 2);
        }
        double yQ_sp = dot(Y.row(i), Q.col(s + p - 1));
        double r = R(s + p - 1, s + p - 1);

        arma::rowvec gammaStar(s + p);

        if(s + p > 1){
          gammaStar.cols(0, s + p - 2) = Y.row(i) * Q.cols(0, s + p - 2);
        }
        gammaStar(s + p - 1) = (pow(yQ_sp, 2) + yNaught2 + r * aNaught * yQ_sp) / (r * aNaught + yQ_sp );

        arma::rowvec solved = solve(R.rows(0, s + p - 1), gammaStar.t()).t();
        BInit(i_vec, pa_i) = solved.cols(s, s + p - 1) * Qp;


        if( s > 0 ){
          OmegaInit(i_vec, sib_i) = solved.cols(0, s - 1);
          OmegaInit(sib_i, i_vec) = solved.cols(0, s - 1).t();
          resid = Y.row(i) - BInit(i_vec, pa_i) * Y.rows(pa_i) - OmegaInit(i_vec, sib_i) * Z.rows(sib2);
        } else {
          resid = Y.row(i) - BInit(i_vec, pa_i) * Y.rows(pa_i);
        }

        arma::vec omegaVec = OmegaInit.col(i);
        omegaVec.shed_row(i);
        OmegaInit(i, i) = 1.0 / Y.n_cols * pow(norm(resid,2), 2) + dot(omegaVec.t() , solve(omegaSubset(i),  omegaVec));
    }

    // succesfully updated. If condition number of Omega[-i,-i] is too large, then it returns 0 and does not update
    return 1;
}


mat lin_sem::omegaSubset(int i)
{
    arma::mat ret;
    ret = OmegaInit;
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

mat lin_sem::getB(){
  return B;
}

mat lin_sem::getOmega(){
  return Omega;
}
