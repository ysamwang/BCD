#ifndef LINSEM_H
#define LINSEM_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp ;
using namespace arma;

class lin_sem
{
public:
    lin_sem(SEXP Br, SEXP Omegar, SEXP BInitr, SEXP OmegaInitr, SEXP Yr, double OmegaInitScale);
    int updateNode(int i, double maxKap);
    void initEst( double OmegaInitScale);
    mat getSigma();
    mat getBInit();
    mat getOmegaInit();
    int getV();
    int singleUpdateOnly(int i);

protected:
    mat Y; // p by n data matrix
    mat Omega; // matrix of 1's and 0's representing bi-directed edge structure
    mat B; // matrix of 1's and 0's representing directed edge structure
    mat BInit; // matrix of directed egdge weights
    mat OmegaInit; //matrix of bidirected edge weights
    mat omegaSubset(int i); // helper function to return Omega[-i,-i]
    mat eyeBSubset(int i); // helper function to return (I - B)[-i, ] 
    int V; // number of nodes 
    // vector which records which nodes only need to be updated once.
    // a 1 in the ith spot indicates the ith node only needs to be updated once
    // this is filled in when updateNode(i) is called
    vec singleUpdates; 
    

private:
};
#endif // LINSEM_H
