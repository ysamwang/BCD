#include "BCD.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List bcdC(SEXP Br, SEXP Omegar, SEXP BInitr, SEXP OmegaInitr, SEXP Yr, int maxIter, int sigConv,
                double maxKap, double tol) {

    //initialize graph object
    lin_sem graph = lin_sem(Br, Omegar, BInitr, OmegaInitr, Yr);
    //update
    updateEstimates(graph, maxIter, sigConv, tol);
    return graph.returnGraph();
}


int updateEstimates(lin_sem graph, int maxIter, int sigConv, double tol) {
    int i;
    int counter = 0;
    double convCrit = 1.0;
    if(sigConv){
        arma::mat oldSigma;
        arma::mat newSigma = graph.getSigma();
        while(convCrit > tol && counter < maxIter){
            oldSigma = newSigma;
            for(i = 0; i < graph.getV(); i++)
            {
                graph.updateNode(i);
            }
            newSigma = graph.getSigma();
            convCrit = norm(oldSigma - newSigma, 2);
            counter ++;
        }
    } else {
        arma::mat oldB;
        arma::mat newB = graph.getBInit();
        arma::mat oldOmega;
        arma::mat newOmega = graph.getOmegaInit();
        while(convCrit > tol && counter < maxIter){
            oldB = newB;
            oldOmega = newOmega;
            for(i = 0; i < graph.getV(); i++)
            {
                graph.updateNode(i);
            }
            newB = graph.getBInit();
            newOmega = graph.getOmegaInit();
            convCrit = norm(oldB - newB, 2) + norm(oldOmega - newOmega, 2);
            counter ++;
        }

    }
  if(counter == maxIter) {
    return 0;
  } else {
    return 1;
  }
}
