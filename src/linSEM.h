#ifndef LINSEM_H
#define LINSEM_H


#include <stdlib.h>
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp ;
using namespace arma;


// lin_sem object holds a mixed graph object for a linear structural equation model
// consisting of directed and bidirected edges. See DFW (INSERT YR)
class lin_sem
{
public:
    // Initialize lin_sem object:
    //
    // 
    //
    // Br is a p x p matrix of 1/0. B_ij = 1 indicates a directed edge j -> i
    // Omegar is a p x p  matrix of 1/0. Omega_ij indicates a bidirected edge i <->j.
    //     Note that Omega_ii must = 1 for all i
    // BInitr is a p x p matrix which holds initial edge weights for the directed edges
    //     If BInitr is null, then the default initialization routine is used
    // OmegaInitr is a p x p matrix which holds initial edge weights for the bidirected edges
    // Yr is a p x n matrix with the observed data. Each row should have mean 0
    // OmegaInitScale is the scaling factor used in the initialization if necessary. See initEst
    lin_sem(SEXP Br, SEXP Omegar, SEXP BInitr, SEXP OmegaInitr, SEXP Yr, double OmegaInitScale);
  
    // update node i through BCD algorithm
    // 
    // i node to be update
    // maxKap is the maximum condition number of OmegaInit[-i,-i]. If this is exceeded,
    //     the method returns 0 and terminates. 
    int updateNode(int i, double maxKap);
    
    // get initial estimates for BInit and OmegaInit
    // 
    
    // OmegaInitScale
    void initEst( double OmegaInitScale);
    mat getSigma();
    mat getBInit();
    mat getOmegaInit();
    int getV();
    int singleUpdateOnly(int i);
    mat getB();
    mat getOmega();

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
