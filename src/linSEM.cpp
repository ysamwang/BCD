#include "linSEM.h"

using namespace Rcpp ;
using namespace arma ;
lin_sem::lin_sem(SEXP Br, SEXP Omegar, SEXP BInitr, SEXP OmegaInitr, SEXP Yr)
{
    B = as<arma::mat>(Br);
    Omega =  as<arma::mat>(Omegar);
    BInit = as<arma::mat>(clone(BInitr));
    OmegaInit =  as<arma::mat>(clone(OmegaInitr));
    Y = as<arma::mat>(Yr);
    V = B.n_cols;
    P = Y.n_cols;
    if(BInit(0,0) == 1)
    {
        initEst();
    }
}

void lin_sem::initEst()
{
    int i;
    //initialize L and Omega
    uvec i_vec(1);
    mat resid;
    resid.copy_size(Y);
    for(i = 0; i < V; i++){

        // Set L through OLS
        uvec pa_i = find(B.row(i));
        vec coef = solve(Y.rows(pa_i), Y.row(i) );
        i_vec[0] = i;
        BInit(i_vec, pa_i ) = coef.t();

        //Estimate residual from OLS and use as omega.ii
        resid.row(i) = (Y.row(i) - coef * Y.rows(pa_i));

        //Randomly initialize remaining elements from std normal
        //Normalize and ensure sum is smaller than sigma 2 to remain diagonally dominant
        OmegaInit = Omega % cov(resid.t());
    }
}


void lin_sem::updateNode(int i)
{
    int k;
    double det1, det2, aNaught;
    // get parents and siblings
    uvec pa_i = find(B.row(i));
    uvec sib_i = find(Omega.row(i));
    uvec i_vec; 
    i_vec<< i <<endr;
    int s = sib_i.n_elem;
    int p = pa_i.n_elem;

    //calculate a_0 and a_pa
    arma::vec a_pa(pa_i.n_elem);

    for(k = 0; k < pa_i.n_elem; k ++){
        //potential speed up
        BInit(i, pa_i(k)) = 2;
        det2 = det(eye(V,V) - B);
        BInit(i, pa_i(k)) = 2;
        det1 = det(eye(V,V) - B);
        a_pa(k) = det2 - det1;
    }
    BInit.row(i) = vec(V, fill::zeros).t();
    aNaught = det(eye(V, V) - B );
    if(norm(a_pa,2) == 0){

    } else {
    //get Qp
    a_pa = a_pa / norm(a_pa);
    aNaught = aNaught / norm(a_pa);

    arma::mat Qp =  join_vert(null(a_pa.t()), a_pa.t());

    // Create pseudo-residuals
    arma::mat Z = solve(omegaSubsetZ(i), eyeBSubset(i));

    // combine into a single matrix
    arma::mat X = join_vert(Z, Qp * Y.rows(pa_i));

    // QR decomposition of X^t
    arma::mat Q, R;
    qr(Q, R, X.t());

    // Compute y_0^2 for gamma^\star
    double yNaught2 = 0;
    for(k = s + p; k < Y.n_cols; k++){
        yNaught2 += pow(dot(Y.row(i), Q.col(k)), 2);
    }
    
    double yQ_sp = dot(Y.row(i), Q.col(s + p));
    double r = R(s+p, s + p);
    arma::rowvec gammaStar(s + p);
    gammaStar.cols(0, s + p - 2) = Y.row(i) * Q.cols(0, s + p - 2);
    gammaStar(s + p) = (pow(yQ_sp,2) + yNaught2 + r * aNaught * yQ_sp) / (r * aNaught + yQ_sp );

    arma::rowvec solved = solve(R.rows(0, s + p - 1), gammaStar.t()).t();

    OmegaInit(i_vec, sib_i) = solved.cols(0, s - 1);
    OmegaInit(sib_i, i_vec) = solved.cols(0, s - 1).t();

    solved = solved * Q;

    BInit(i_vec, pa_i) = solved.cols(s, s + p - 1);
    arma::rowvec resid = Y.row(i) - BInit(i_vec, pa_i) * Y(i_vec, pa_i) -  OmegaInit(i_vec, sib_i) * Z;

    arma::vec omegaVec = Omega.col(i);
    omegaVec.shed_row(i);
    OmegaInit(i, i) = 1.0 / Y.n_cols * pow (norm(resid), 2) +
         dot(omegaVec.t() , solve(omegaSubset(i),  omegaVec));

    }
}

mat lin_sem::omegaSubsetZ(int i)
{
    arma::mat ret;
    if(i == 0){
        ret =  OmegaInit(span(1,V), span(1,V));
    } else if (i == V) {
        ret = OmegaInit(span(0, V-1), span(0, V-1));
    } else {
        ret = OmegaInit;
        ret.shed_row(i);
        ret.shed_col(i);
    }
    arma::uvec zsib_ind = find(Omega.row(i));
    zsib_ind(find(zsib_ind > i )) = zsib_ind(find(zsib_ind > i )) - 1;
    return ret.rows(zsib_ind);
}

mat lin_sem::omegaSubset(int i)
{
    arma::mat ret;
    if(i == 0){
        ret =  OmegaInit(span(1,V - 1), span(1, V -1));
    } else if (i == V) {
        ret = OmegaInit(span(0, V - 2), span(0, V -2));
    } else {
        ret = OmegaInit;
        ret.shed_row(i);
        ret.shed_col(i);
    }
    return ret;
}


mat lin_sem::eyeBSubset(int i)
{
    arma::mat ret = eye(V, V) - B;
    ret.shed_row(i);
    return ret;
}


mat lin_sem::getBInit(){
    return BInit;
}

mat lin_sem::getOmegaInit(){
    return OmegaInit;
}

mat lin_sem::getSigma() {
    return solve(eye(V, V) - B, solve(eye(V, V) - BInit, OmegaInit).t()).t();
}

Rcpp::List lin_sem::returnGraph() {
        return Rcpp::List::create(Rcpp::Named("OmegaInit", OmegaInit),
                              Rcpp::Named("BInit", BInit));
}

int lin_sem::getV()
{
    return V;
}

