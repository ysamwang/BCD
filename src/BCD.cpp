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
          // check if node needs to be updated
          // after the first pass, if there are no incident bidirected edges
          // and no parents in a directed cycle, then the node does not need to be
          // updated on subsequent passes and graph.SingleUpdateOnly(i) will be set to 1
          if(!graph.singleUpdateOnly(i)){

            // Update node returns a 0 if the update is unsuccesful
            // ie when the conditioning number is too large
            // in this case, simply return current estimates
            if(!graph.updateNode(i, maxKap)){
              return Rcpp::List::create(Rcpp::Named("SigmaHat", graph.getSigma()),
                                      Rcpp::Named("OmegaHat", graph.getOmegaInit()),
                                      Rcpp::Named("BHat", graph.getBInit()),
                                      Rcpp::Named("Iter", counter),
                                      Rcpp::Named("Converged", 0));
          }
        }
      }

      //update convergence criteria and counter
      convCrit =  accu(abs(oldSigma - graph.getSigma())) / accu(abs(oldSigma));
      counter ++;
    }
  } else {
    //place holders for old B and old Omega
    arma::mat oldB;
    arma::mat oldOmega;
    while(convCrit > tol && counter < maxIter){
      oldB = graph.getBInit();
      oldOmega = graph.getOmegaInit();


      //update nodes
      for(i = 0; i < graph.getV(); i++)
      {
        // check if node needs to be updated
        // after the first pass, if there are no incident bidirected edges
        // and no parents in a directed cycle, then the node does not need to be
        // updated on subsequent passes and graph.SingleUpdateOnly(i) will be set to 1
        if(!graph.singleUpdateOnly(i)){

          // Update node returns a 0 if the update is unsuccesful
          // ie when the conditioning number is too large
          // in this case, simply return current estimates
          if(!graph.updateNode(i, maxKap)){
            return Rcpp::List::create(Rcpp::Named("SigmaHat", graph.getSigma()),
                                      Rcpp::Named("OmegaHat", graph.getOmegaInit()),
                                      Rcpp::Named("BHat", graph.getBInit()),
                                      Rcpp::Named("Iter", counter),
                                      Rcpp::Named("Converged", 0));
          }
        }
      }

    convCrit =  (accu(abs(oldB - graph.getBInit())) + accu(abs(oldOmega - graph.getOmegaInit()))) / (accu(graph.getB()) + accu(graph.getOmega()));
    counter++;
    } //end while

  }

    //while statement terminates, return current estimates
    return Rcpp::List::create(Rcpp::Named("SigmaHat", graph.getSigma()),
                              Rcpp::Named("OmegaHat", graph.getOmegaInit()),
                              Rcpp::Named("BHat", graph.getBInit()),
                              Rcpp::Named("Iter", counter),
                              Rcpp::Named("Converged", (convCrit < tol)));
}
