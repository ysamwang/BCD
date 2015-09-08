#include "BCD.h"

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List bcdC(SEXP Br, SEXP Omegar, SEXP BInitr, SEXP OmegaInitr, SEXP Yr, int maxIter, int sigConv,
                double maxKap, double tol, double omegaInitScale) {

    //initialize graph object
  lin_sem graph = lin_sem(Br, Omegar, BInitr, OmegaInitr, Yr, omegaInitScale);
  
  //start Updates
  int i;
  int counter = 0;
  double convCrit = 1.0;

  if(sigConv){
    
    arma::mat oldSigma; // hold previous sigma estimate

    while(convCrit > tol && counter < maxIter){
      oldSigma = graph.getSigma();
      
      //update nodes
      for(i = 0; i < graph.getV(); i++)
      {
          if(!graph.updateNode(i, maxKap)){
            return Rcpp::List::create(Rcpp::Named("SigmaHat", graph.getSigma()),
                                      Rcpp::Named("OmegaHat", graph.getOmegaInit()),
                                      Rcpp::Named("BHat", graph.getBInit()),
                                      Rcpp::Named("Iter", counter),
                                      Rcpp::Named("Converged", 0)
            );
          }
      }
      
      //update convergence criteria and counter
      convCrit = norm(oldSigma - graph.getSigma(), "fro");
      counter ++;
    }
  } else {
    //place holders for old B and old Omega
    arma::mat oldB;
    arma::mat oldOmega;
    while(convCrit > tol && counter < maxIter){
      oldB = graph.getBInit();
      oldOmega = graph.getOmegaInit();
      for(i = 0; i < graph.getV(); i++)
      {
        //check if the node needs to be updated
          if(!graph.updateNode(i, maxKap)){
            return Rcpp::List::create(Rcpp::Named("SigmaHat", graph.getSigma()),
                                    Rcpp::Named("OmegaHat", graph.getOmegaInit()),
                                    Rcpp::Named("BHat", graph.getBInit()),
                                    Rcpp::Named("Iter", counter),
                                    Rcpp::Named("Converged", 0)
          );
        }
      }
     
     convCrit = norm(oldB - graph.getBInit(), "fro") + norm(oldOmega - graph.getOmegaInit(), "fro");
    counter++;
    }
  }

    return Rcpp::List::create(Rcpp::Named("SigmaHat", graph.getSigma()),
                              Rcpp::Named("OmegaHat", graph.getOmegaInit()),
                              Rcpp::Named("BHat", graph.getBInit()),
                              Rcpp::Named("Iter", counter),
                              Rcpp::Named("Converged", (convCrit < tol))
                                );
}
