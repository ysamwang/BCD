#' Fitting Linear Structural Equation Models
#' 
#' 
#' Fits Empirical Likelihood estimates for given SEM
#' 
#' @param Y p x n matrix with data
#' @param B V by V matrix with {0,1} giving structure of directed edges
#' @param covarRestrict
#' @param meanEst
#' @param tol
#' @param outerTol
sempl_cd <- function(Y, B, Omega, covarRestrict, meanEst = F, tol = 1e-3, outerTol = 1e-2){

  # number of nodes
  V <- dim(B)[1]
  n <- dim(Y)[2]

  ## If fitting a completely empty model
  if(meanEst == F & sum(B) == 0){
    
    mod <- sempl_input_weights(b_weights_r = c(), y_r = Y, b_r = B, mean_est_r = meanEst, covar_restrict = covarRestrict,
                               tol = tol, max_iter = 100)
    return(sum(log(mod$d)) - n * log(n))
    
  } else { 
    
    out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y - rowMeans(Y), BInit = NULL,
                            OmegaInit = NULL, sigConv = 0, maxIter = 0,
                            msgs = FALSE, omegaInitScale = .9)
    
    B.hat <- out_ricf_init$BHat
  
    ## Set mu.hat if necessary
    if (meanEst) { 
      mu.hat <- solve(diag(rep(1,V)) - B.hat, rowMeans(Y)) 
    } else {
      mu.hat <- NULL
    }
    
    
    ## If mu.hat is not needed, it should still be null
    init_val <- c(mu.hat, B.hat[which(B == 1)])
    
    out <- nlm(f = sempl_input, p = init_val, y_r = Y, b_r = B,
               mean_est_r = meanEst, covar_restrict = covarRestrict,
               tol = tol, max_iter = 100, check.analyticals = F,
               gradtol = outerTol)
    return(out$minimum - n * log(n))
  }
}

