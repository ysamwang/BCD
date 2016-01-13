fitEL <- function(Y, B, Omega, B.hat = NULL, Omega.hat = NULL, mu.hat = NULL, method = "BFGS",
                  tol = 1e-6, maxInnerIter = 100, maxOuterIter = 100, meanEst = 0, naive = FALSE, high_moments = 0){
  V <- dim(B)[1]
  
  ### Initialization points ###
  
  # For B
  
  
  if(!naive){
    if(is.null(B.hat)){
      out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y - rowMeans(Y), BInit = NULL,
                            OmegaInit = NULL, sigConv = 0, maxIter = 0,
                            msgs = FALSE, omegaInitScale = .9)
      init_val = c(out_ricf_init$BHat[B==1])
    } else {
      init_val = c(B.hat[B==1])
    }
  # For the mean (if needed)0
  if(meanEst & !is.null(mu.hat)){
    init_val <- c(mu.hat, init_val)
  } else if(meanEst){
    init_val <- c( (diag(rep(1,V)) - out_ricf_init$BHat) %*% rowMeans(Y), init_val)
  }
  
  # Number of restrictions
    if(high_moments){
      num_dual_vars = V + 3 * sum(Omega == 0)/2
    } else {
      num_dual_vars = V + sum(Omega == 0)/2
    }
  
  
  
    optim_out <- optim(par = init_val, fn = sem_el_fit_obj, gr = NULL,
                       y_r = Y, omega_r = Omega, b_r = B,
                       dual_r = rep(0, num_dual_vars), tol = tol, max_iter = maxInnerIter,
                       meanEst = meanEst, high_moments = high_moments, method = method, control = list(fnscale = -1))

    fitted_mod <- sem_el_fit_weights(optim_out$par, y_r = Y, omega_r = Omega, b_r = B,
                                   dual_r = rep(0, num_dual_vars), tol = tol,
                                   max_iter = maxInnerIter, meanEst = meanEst, high_moments = high_moments)
    B.hat <- matrix(0, nrow = V, ncol = V)
    
    p <- c(1/fitted_mod$d)
    
    if(meanEst){
      mu.hat <- optim_out$par[c(1:V)]
      B.hat[B==1] <- optim_out$par[-c(1:V)]
      temp <- (diag(rep(1,V)) - B.hat) %*% (Y - solve(diag(rep(1,V)) - B.hat, mu.hat))  %*% diag(sqrt(p))
      
    } else {
      B.hat[B==1] <- optim_out$par
      temp <- (diag(rep(1,V)) - B.hat) %*% Y  %*% diag(sqrt(p))
    }

    
    

    Omega.hat <-  round(temp %*% t(temp),10)
    lrt <- -2 * sum(log(p * n))
    
    
  } else {  #### Naive Method ####
   
    
    ## Initialize guesses
    if(is.null(B.hat)){
      out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                            OmegaInit = NULL, sigConv = 0, maxIter = 0,
                            msgs = FALSE, omegaInitScale = .9)
      init_val = c(out_ricf_init$BHat[B==1], out_ricf_init$OmegaHat[lower.tri(Omega, diag = T) & Omega == 1])
    } else {
      init_val = c(B.hat[B==1], Omega.hat[Omega == 1])
    }
    # For the mean (if needed)
    if(meanEst & !is.null(mu.hat)){
      init_val <- c(mu.hat, init_val)
    } else if(meanEst){
      init_val <- c(rowMeans(Y), init_val)
    }
    
    ## Fit Model ##
    num_dual_vars = V + V * (V + 1)/2
    
    
    optim_out <- optim(par = init_val, fn = sem_el_naive_fit_obj, gr = NULL,
                       y_r = Y, omega_r = Omega, b_r = B,
                       dual_r = rep(0, num_dual_vars), tol = tol, max_iter = maxInnerIter,
                       meanEst = meanEst, method = method, control = list(fnscale = -1))
    
    fitted_mod <- sem_el_naive_fit_weights(optim_out$par,y_r = Y, omega_r = Omega, b_r = B,
                                     dual_r = rep(0, num_dual_vars), tol = tol, max_iter = maxInnerIter, meanEst = meanEst)
    
    B.hat <- matrix(0, nrow = V, ncol = V)
    
    if(meanEst){
      mu.hat <- optim_out$par[c(1:V)]
      B.hat[B==1] <- optim_out$par[(c(1:sum(B)) + V)]
      Omega.hat[Omega == 1 & lower.tri(Omega, diag = T)] <- optim_out$par[-c(1:(V + sum(B)))]
      
    } else {
      B.hat[B==1] <- optim_out$par[c(1:sum(B))]
      Omega.hat[Omega == 1 & lower.tri(Omega, diag = T)] <- optim_out$par[-c(1:sum(B))]
    }
    Omega.hat <- Omega.hat + t(Omega.hat) - diag(diag(Omega.hat))
    
    p <- c(1/fitted_mod$d)
    lrt <- -2 * sum(log(p * n))
    
  } # end naive

  return(list(B.hat = B.hat,
              Omega.hat = Omega.hat,
              mu.hat = mu.hat,
              p = p,
              converged = !optim_out$convergence & (abs(sum(p)-1) < 1e-8),
              lrt = lrt)
         )
}



# Calculate weights at a specified point
fitEL.weights <- function(Y, B, Omega, B.hat, mu.hat = NULL,
                  tol = 1e-6, maxInnerIter = 100){

  num_dual_vars = V + sum(Omega == 0)/2
  
  val <- B.hat[which(B==1)]
  
  fitted_mod <- sem_el_fit_weights(val, y_r = Y, omega_r = Omega, b_r = B,
                                   dual_r = rep(0, num_dual_vars), tol = tol, max_iter = maxInnerIter, meanEst = !is.null(mu.hat))
  
  p <- c(1/fitted_mod$d)
  temp <- (diag(rep(1,V)) - B.hat) %*% Y %*% diag(sqrt(p))
  Omega.hat <-  round(temp %*% t(temp),10)
  lrt <- -2 * sum(log(p * n))

  return(list(Omega.hat = Omega.hat,
              p = p,
              converged = (abs(sum(p)-1) < 1e-8),
              lrt = lrt))
  
}