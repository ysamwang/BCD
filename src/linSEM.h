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
    //constructor
    lin_sem(SEXP Br, SEXP Omegar, SEXP BInitr, SEXP OmegaInitr, SEXP Yr);
    void updateNode(int i);
    void initEst();
    mat getSigma();
    mat getBInit();
    mat getOmegaInit();
    int getV();
    int getP();
    Rcpp::List returnGraph();

protected:

    mat Y;
    mat Omega;
    mat B;
    mat BInit;
    mat OmegaInit;
    mat omegaSubsetZ(int i);
    mat omegaSubset(int i);
    mat eyeBSubset(int i);
    int V;
    int P;

private:
};
#endif // LINSEM_H
