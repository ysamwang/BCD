#' Fitting Linear Structural Equation Models
#' 
#' 
#' Fits Empirical Likelihood estimates for given SEM
#' 
#' @param Y p x n matrix with data
#' @param B V by V matrix with {0,1} giving structure of directed edges
#' @param Omega V by V matrix with {0,1} giving structure of bi-directed edges
#' @param B.hat V by V matrix giving initial edges weights for directed edges. If BInit is NULL,
#'    a default initialization will be used. 
#' @param Omega.hat V by V matrix giving initial edge weights for bi-directed edges
#' @param mu.hat Initial value for means  
#' @param tol
#' @param maxInnerIter
#' @param outerTol
#' @param meanEst
#' @param type profile, naive, euclid
#' @param covarRestrict
#' @export
sempl <- function(Y, B, Omega, B.hat = NULL, Omega.hat = NULL, mu.hat = NULL, meanEst = T,
                  covarRestrict = NULL, type = "profile", tol = 1e-4, maxIter = 500,
                  outerIters = 200){
  
  ##############
  # Initialize #
  ##############
  
  
  if(type == "profile"){
    
    if(is.null(covarRestrict)){
      temp <- which(Omega == 0, arr.ind = T)
      covarRestrict <- unique(cbind(apply(temp, MAR = 1, max), apply(temp, MAR = 1, min))) - 1
    }
    
    # number of nodes
    V <- dim(B)[1]
    n <- dim(Y)[2]
    
    ## If fitting a completely empty model
    if(meanEst == F & sum(B) == 0){
      mod <- sempl_input_weights(b_weights_r = c(), y_r = Y, b_r = B, mean_est_r = meanEst, covar_restrict = covarRestrict,
                                 tol = tol, max_iter = maxIter)
      p <- c(1/mod$d)
      return(list(mu.hat = NULL,
                  B.hat = matrix(0, nrow = V, ncol = V),
                  Omega.hat = Omega.hat,
                  p = p,
                  lr =  mod$lr - n * log(n),
                  converged = abs(sum(p) - 1) <1e-6,
                  iter = 1))
    } else { ## If fitting non-empty model
    
      ## Using RICF initialization if necessary
      if(is.null(B.hat) ) {
        out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y - rowMeans(Y), BInit = NULL,
                              OmegaInit = NULL, sigConv = 0, maxIter = 1,
                              msgs = FALSE, omegaInitScale = .9)
        
      }
      
      ## Set B.hat if necessary
      if(is.null(B.hat)){
        B.hat <- out_ricf_init$BHat
      }
      
      ## Set mu.hat if necessary
      if((is.null(mu.hat) & meanEst)){
        mu.hat <- solve(diag(rep(1,V)) - B.hat, rowMeans(Y))
      }
      
      
      ## If mu.hat is not needed, it should still be null
      init_val <- c(mu.hat, B.hat[which(B == 1)])
      
      out <- nlm(f = sempl_input, p = init_val, y_r = Y, b_r = B,
                 mean_est_r = meanEst, covar_restrict = covarRestrict,
                 tol = tol, max_iter = maxIter, check.analyticals = F,
                 iterlim = outerIters)
      
      mod <- sempl_input_weights(out$estimate, Y, B, meanEst, covarRestrict,
                                 tol, maxIter)
      
      p <- c(1/mod$d)
      ## Get fitted values
      if(meanEst){
        mu.hat <- out$estimate[1:V]
        B.hat <- matrix(0, nrow = V, ncol = V)
        B.hat[which(B==1)] <- out$estimate[-c(1:V)]
        temp <- (diag(rep(1,V)) - B.hat) %*% Y
        Omega.hat <- temp %*% diag(p) %*% t(temp)
      } else {
        B.hat <- matrix(0, nrow = V, ncol = V)
        B.hat[which(B==1)] <- out$estimate
        temp <- (diag(rep(1,V)) - B.hat) %*% Y
        Omega.hat <- (temp %*% diag(p) %*% t(temp)) * Omega
      }
      
      return(list(mu.hat = mu.hat,
                  B.hat = B.hat,
                  Omega.hat = Omega.hat,
                  p = p,
                  lr =  mod$lr - n * log(n),
                  converged = out$code,
                  iter = out$iter))
    }
    
    
  } else if(type == "naive"){
    
    # number of nodes
    V <- dim(B)[1]
    n <- dim(Y)[2]
    
    ## Using RICF initialization if necessary
    if( is.null(B.hat) | is.null(Omega.hat) ) {
      out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y - rowMeans(Y), BInit = NULL,
                            OmegaInit = NULL, sigConv = 0, maxIter = 0,
                            msgs = FALSE, omegaInitScale = .9)
      
    }
    
    ## Set B.hat if necessary
    if(is.null(B.hat)){
      B.hat <- out_ricf_init$BHat
    }
    
    ## Set Omega.hat if necessary
    if( is.null(Omega.hat) ){ 
      Omega.hat <- out_ricf_init$OmegaHat
    }
    
    ## Set mu.hat if necessary
    if((is.null(mu.hat) & meanEst)){
      mu.hat <- solve(diag(rep(1,V)) - B.hat, rowMeans(Y))
    }
    
    
    ## If mu.hat is not needed, it should still be null
    init_val <- c(mu.hat, B.hat[which(B == 1)],
                  Omega.hat[(Omega == 1) & lower.tri(Omega, diag = T)])
    
    out <- nlm(f = sempl_input_naive, p = init_val, y_r = Y, b_r = B,
               omega_r = Omega, mean_est_r = meanEst, tol = tol,
               max_iter = maxIter, check.analyticals = F, iterlim = outerIters)
    mod <- sempl_input_naive_weights(out$estimate, Y, B, Omega,
                                     meanEst, tol, maxIter)
    
    p <- c(1/mod$d)
    ## Get fitted values
    if(meanEst){
      mu.hat <- out$estimate[1:V]
      B.hat <- matrix(0, nrow = V, ncol = V)
      B.hat[which(B==1)] <- out$estimate[c(1:sum(B)) + V]
      Omega.hat[(Omega == 1) & lower.tri(Omega, diag = T)] <- out$estimate[-c(1:sum(B) + V)]   
      Omega.hat <- Omega.hat + t(Omega.hat) - diag(diag(Omega.hat))
    } else {
      B.hat <- matrix(0, nrow = V, ncol = V)
      B.hat[which(B==1)] <- out$estimate[c(1:sum(B))]
      Omega.hat <- matrix(0, nrow = V, ncol = V)
      Omega.hat[(Omega == 1) & lower.tri(Omega, diag = T)] <- out$estimate[-c(1:sum(B))]   
      Omega.hat <- Omega.hat + t(Omega.hat) - diag(diag(Omega.hat))
    }
    
    return(list(B.hat = B.hat,
                Omega.hat = Omega.hat,
                p = p,
                lr =  mod$lr - n * log(n),
                converged = out$code,
                iter = out$iter))
  }
  
}