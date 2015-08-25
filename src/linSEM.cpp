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
    //initialize B and Omega
    uvec i_vec(1);
    mat resid;
    resid.copy_size(Y);
    for(i = 0; i < V; i++) {
        // Set B through OLS
        uvec pa_i = find(B.row(i));
        if( pa_i.n_rows > 0) {
            vec coef = solve(Y.rows(pa_i).t(), Y.row(i).t());
            i_vec[0] = i;
            BInit(i_vec, pa_i ) = coef.t();

            //Estimate residual from OLS and use as initial estimates
            resid.row(i) = (Y.row(i) - coef * Y.rows(pa_i));
        } else {
            resid.row(i) = Y.row(i) - mean(Y.row(i));
        }
    }
    OmegaInit = Omega % cov(resid.t());
}


void lin_sem::updateNode(int i)
{
    int k;
    double det1, det2, aNaught;
    // get parents and siblings
    uvec pa_i = find(B.row(i));
    uvec sib_i = find(Omega.row(i));
    arma::uvec i_in_sib_i = find(sib_i == i);
    arma::uvec sib2;
    int i_index = i_in_sib_i[0];
    sib_i.shed_row(i_index);
    
    uvec i_vec;
    i_vec<< i <<endr;
    int s = sib_i.n_elem;
    int p = pa_i.n_elem;
    arma::mat X;
    arma::mat Z;
    arma::vec solved;
    arma::rowvec resid;
    

    //calculate a_0 and a_pa
    arma::vec a_pa(pa_i.n_elem);

    for(k = 0; k < pa_i.n_elem; k ++) {
        //potential speed up
        BInit(i, pa_i(k)) = 2.0;
        det2 = det(eye(V,V) - BInit);
        BInit(i, pa_i(k)) = 1.0;
        det1 = det(eye(V,V) - BInit);
        a_pa[k] = det2 - det1;
    }

    BInit.row(i) = vec(V, fill::zeros).t();
    aNaught = det(eye(V, V) - B );
    
    Rcout << "About to Start BCD: "<< i <<"; Norm is: " <<norm(a_pa,2) <<endl;
    //Acyclic case or case where no parents are involved in cycles
    if(norm(a_pa,2) == 0) {
      if(s > 0 ){
        //get psuedo variables
        arma::uvec sib2 = sib_i;
        sib2.elem(find(sib2 > i)) = sib2.elem(find(sib2 > i)) - 1;
        
        arma::mat Z = solve(omegaSubset(i), eyeBSubset(i)) * Y;
        
        if( p > 0 ) {
        // combine into a single matrix
          X = join_vert(Z.rows(sib2), Y.rows(pa_i));
          solved = solve(X.t(), Y.row(i).t());
          
          //update Binit since there are
          BInit(i_vec, pa_i) = solved.rows(s, s + p - 1).t();
        } else {
          X = Z.rows(sib2);
          solved = solve(X.t(), Y.row(i).t());
        }
        
        //update OmegaInit since there are siblings
        OmegaInit(i_vec, sib_i) = solved.rows(0, s - 1).t();
        OmegaInit(sib_i, i_vec) = solved.rows(0, s - 1);
        
        resid = Y.row(i) - solved.t() * X;
        
      
      } else {
        if( p > 0 ) {
          X = Y.rows(pa_i);
          solved = solve(X.t(), Y.row(i).t());
          BInit(i_vec, pa_i) = solved.cols(s, s + p - 1);
          resid = Y.row(i) - solved.t() * X;
        } else {
          resid = Y.row(i) - mean(Y.row(i));
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
        
        Rcout <<"About to get QP" <<endl;
        a_pa = a_pa / norm(a_pa);
        aNaught = aNaught / norm(a_pa);
        
        Rcout << "About to join" <<endl;
        if(a_pa.n_elem > 1) {
            Rcout <<"join"<<endl <<null(a_pa.t()) << endl << a_pa.t() <<endl;
            Qp =  join_vert(null(a_pa.t()).t(), a_pa.t()); 
        } else {
          Qp << 1.0 <<endr;
        }

        
        if(s > 0){
          //get psuedo variables
          Rcout << "About to get pseudo" <<endl;
          sib2 = sib_i;
          sib2.elem(find(sib2 > i)) = sib2.elem(find(sib2 > i)) - 1;
          Z = solve(omegaSubset(i), eyeBSubset(i)) * Y;
          X = join_vert(Z.rows(sib2), Qp * Y.rows(pa_i));
          
        } else {
          X =  Qp * Y.rows(pa_i);
        }
        Rcout <<"Getting QR"<<endl;
        // QR decomposition of X^t
        arma::mat Q, R;
        qr(Q, R, X.t());
        Rcout << "About to Compute" <<endl;
        // Compute y_0^2 for gamma^\star
        double yNaught2 = 0;
        for(k = s + p; k < Y.n_cols; k++) {
            yNaught2 += pow(dot(Y.row(i), Q.col(k)), 2);
        }
        double yQ_sp = dot(Y.row(i), Q.col(s + p - 1));
        double r = R(s + p - 1, s + p - 1);
        
        arma::rowvec gammaStar(s + p);

        Rcout << "About to multiply" <<endl;
        if(s + p > 1){ 
          gammaStar.cols(0, s + p - 2) = Y.row(i) * Q.cols(0, s + p - 2);
        }
        
        gammaStar(s + p - 1) = (pow(yQ_sp, 2) + yNaught2 + r * aNaught * yQ_sp) / (r * aNaught + yQ_sp );
        Rcout << "final step" <<endl;
        arma::rowvec solved = solve(R.rows(0, s + p - 1), gammaStar.t()).t();
        
        Rcout << "About to update B" <<endl;
        BInit(i_vec, pa_i) = solved.cols(s, s + p - 1) * Qp;
        
        if( s > 0 ){
          Rcout << "About to Update Omega" << endl;
          OmegaInit(i_vec, sib_i) = solved.cols(0, s - 1);
          OmegaInit(sib_i, i_vec) = solved.cols(0, s - 1).t();
          resid = Y.row(i) - BInit(i_vec, pa_i) * Y.rows(pa_i) - OmegaInit(i_vec, sib_i) * Z.rows(sib2);
        } else {
          resid = Y.row(i) - BInit(i_vec, pa_i) * Y.rows(pa_i);
        }

        
        Rcout <<"Getting var Est" <<endl;
        arma::vec omegaVec = OmegaInit.col(i);
        omegaVec.shed_row(i);
        OmegaInit(i, i) = 1.0 / Y.n_cols * pow (norm(resid), 2) +
                          dot(omegaVec.t() , solve(omegaSubset(i),  omegaVec));

        Rcout << "Estimated Var: "<< OmegaInit(i, i) <<std::endl;
    }
}


mat lin_sem::omegaSubset(int i)
{
    arma::mat ret;
    if(i == 0) {
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
    return solve(eye(V, V) - B, solve(eye(V, V) - BInit, OmegaInit).t()).t();
}

int lin_sem::getV()
{
    return V;
}

