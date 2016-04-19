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
#' @param overIdRestrict
#' @export
sempl <- function(Y, B, Omega, B.hat = NULL, Omega.hat = NULL, mu.hat = NULL, meanEst = T,
                  overIdRestrict = NULL, type = "profile"){
  
  ##############
  # Initialize #
  ##############
  
  # number of nodes
  V <- dim(B)[1]
  
  ## Using RICF initialization if necessary
  if(is.null(B.hat) | (is.null(Omega.hat) & type == "naive")) {
      out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y - rowMeans(Y), BInit = NULL,
                            OmegaInit = NULL, sigConv = 0, maxIter = 0,
                            msgs = FALSE, omegaInitScale = .9)

  }
  
  ## Set B.hat if necessary
  if(is.null(B.hat)){
    B.hat <- out_ricf_init$bHat
  }
  
  ## Set Omega.hat if necessary
  if((is.null(Omega.hat) & type == "naive")){
    Omega.hat <- out_ricf_init$omegaHat
  }
  
  ## Set mu.hat if necessary
  if((is.null(mu.hat) & meanEst)){
    mu.hat <- solve(diag(rep(1,V)) - B.hat, rowMeans(Y))
  }
  
  


  if(type == "profile"){
    
    
    
    if(is.null(B.hat)){
        out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y - mu.hat, BInit = NULL,
                            OmegaInit = NULL, sigConv = 0, maxIter = 0,
                            msgs = FALSE, omegaInitScale = .9)
      init_val <- c(out_ricf_init$BHat[B==1])
    } else {
      init_val <- c(B.hat[B == 1])
    }
    
    
    # For the mean (if needed)
    if(meanEst){
      init_val <- c(solve(rep(1, V) - B.hat) %*% mu.hat, init_val)
    }
    
    moment_2_restrictions_r <- which(Omega == 0 & lower.tri(Omega), arr.ind = TRUE) - 1
    if(high_moments == 3){
      
      moment_3_restrictions_r <- rbind(cbind(moment_2_restrictions_r, moment_2_restrictions_r[,1]),
                                       cbind(moment_2_restrictions_r, moment_2_restrictions_r[,2]))
      moment_4_restrictions_r <- NULL
    } else if (high_moments == 4) {
      moment_3_restrictions_r <- rbind(cbind(moment_2_restrictions_r, moment_2_restrictions_r[,1]),
                                       cbind(moment_2_restrictions_r, moment_2_restrictions_r[,2]))
      
      moment_4_restrictions_r <- rbind(cbind(moment_2_restrictions_r, moment_2_restrictions_r[,1], moment_2_restrictions_r[,1]),
                                       cbind(moment_2_restrictions_r, moment_2_restrictions_r[,2],  moment_2_restrictions_r[,2]))
    } else {
      moment_3_restrictions_r <- NULL
      moment_4_restrictions_r <- NULL
    }
    
    
    optim_out <- nlm(f = sem_el_fit_obj_nlm, p = init_val,
                       y_r = Y, b_r = B, moment_2_restrictions_r = moment_2_restrictions_r,
                       moment_3_restrictions_r = moment_3_restrictions_r, moment_4_restrictions_r = moment_4_restrictions_r,
                       tol = tol, max_iter = maxInnerIter, meanEst = meanEst,)
    
    fitted_mod <- sem_el_fit_weights(optim_out$par, y_r = Y, b_r = B, moment_2_restrictions_r = moment_2_restrictions_r,
                                     moment_3_restrictions_r = moment_3_restrictions_r, moment_4_restrictions_r = moment_4_restrictions_r,
                                     tol = tol, max_iter = maxInnerIter, meanEst = meanEst)
    B.hat <- matrix(0, nrow = V, ncol = V)
    
    p <- c(1/fitted_mod$d)
    
    if(meanEst){
      mu.hat <- optim_out$par[c(1:V)]
      B.hat[B==1] <- optim_out$par[-c(1:V)]
      temp <- (diag(rep(1,V)) - B.hat) %*% (Y - solve(diag(rep(1,V)) - B.hat, mu.hat))  %*% diag(sqrt(p))
      
    } else {
      B.hat[B == 1] <- optim_out$par
      temp <- (diag(rep(1,V)) - B.hat) %*% Y  %*% diag(sqrt(p))
    }
    
    
    
    
    Omega.hat <-  round(temp %*% t(temp),10)
    lrt <- -2 * sum(log(p * n))
    
    temp <- solve(diag(rep(1,V)) - B.hat)
    
    return(list(B.hat = B.hat,
                Omega.hat = Omega.hat,
                mu.hat = mu.hat,
                p = p,
                dual = fitted_mod$dual,
                converged = !optim_out$convergence & (abs(sum(p)-1) < 1e-8),
                lrt = lrt,
                Sigma.hat = temp %*% Omega.hat %*% t(temp)))
    
    
  } else if (type == "naive"){  #### Naive Method ####
    V <- dim(B)[1]
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
    
    optim_out <- optim(par = init_val, fn = sem_el_naive_fit_obj, gr = NULL,
                       y_r = Y, omega_r = Omega, b_r = B, dual_r = rep(0, V + V * (V+1)/2), tol = tol,
                       max_iter = maxInnerIter, meanEst = meanEst, method = optim.method, control = list(fnscale = -1, reltol = outerTol))
    
    
    
    B.hat <- Omega.hat <- matrix(0, nrow = V, ncol = V)
    
    if(meanEst){
      mu.hat <- optim_out$par[c(1:V)]
      B.hat[B==1] <- optim_out$par[(c(1:sum(B)) + V)]
      Omega.hat[Omega == 1 & lower.tri(Omega, diag = T)] <- optim_out$par[-c(1:(V + sum(B)))]
      
    } else {
      B.hat[B==1] <- optim_out$par[c(1:sum(B))]
      Omega.hat[Omega == 1 & lower.tri(Omega, diag = T)] <- optim_out$par[-c(1:sum(B))]
    }
    Omega.hat <- Omega.hat + t(Omega.hat) - diag(diag(Omega.hat))
    
    
    fitted_mod <- sem_el_naive_fit_weights(optim_out$par, y_r = Y, omega_r = Omega, b_r = B,
                                           dual_r = rep(0, V + V * (V+1)/2), meanEst = meanEst, tol = tol,
                                           max_iter = maxInnerIter) 
    
    
    p <- c(1/fitted_mod$d)
    lrt <- -2 * sum(log(p * n))
    temp <- solve(diag(rep(1,V)) - B.hat)

        return(list(B.hat = B.hat,
                Omega.hat = Omega.hat,
                mu.hat = mu.hat,
                p = p,
                converged = !optim_out$convergence & (abs(sum(p)-1) < 1e-8),
                lrt = lrt,
                Sigma.hat = temp %*% Omega.hat %*% t(temp)))
    
  } else if(type == "euclid"){
    
    V <- dim(B)[1]
    ## Initialize guesses
    if(is.null(B.hat)){
      if(meanEst){
      out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y - rowMeans(Y), BInit = NULL,
                            OmegaInit = NULL, sigConv = 0, maxIter = 0,
                            msgs = FALSE, omegaInitScale = .9)
      } else {
        out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                              OmegaInit = NULL, sigConv = 0, maxIter = 0,
                              msgs = FALSE, omegaInitScale = .9)
      }
      init_val = c(out_ricf_init$BHat[B==1])
    } else {
      init_val = c(B.hat[B == 1])
    }
    
    # For the mean (if needed)
    if(meanEst & !is.null(mu.hat)){
      init_val <- c(mu.hat, init_val)
    } else if(meanEst){
      if(is.null(B.hat)){
        init_val <- c((diag(rep(1,V)) - out_ricf_init$BHat) %*% rowMeans(Y), init_val)
      } else {
        init_val <- c((diag(rep(1,V)) - B.hat) %*% rowMeans(Y), init_val)
      }
    }
    
    
    moment_2_restrictions_r <- which(Omega == 0 & lower.tri(Omega), arr.ind = TRUE) - 1
    if(high_moments == 3){
      
      moment_3_restrictions_r <- rbind(cbind(moment_2_restrictions_r, moment_2_restrictions_r[,1]),
                                       cbind(moment_2_restrictions_r, moment_2_restrictions_r[,2]))
      moment_4_restrictions_r <- NULL
    } else if (high_moments == 4) {
      moment_3_restrictions_r <- rbind(cbind(moment_2_restrictions_r, moment_2_restrictions_r[,1]),
                                       cbind(moment_2_restrictions_r, moment_2_restrictions_r[,2]))
      
      moment_4_restrictions_r <- rbind(cbind(moment_2_restrictions_r, moment_2_restrictions_r[,1], moment_2_restrictions_r[,1]),
                                       cbind(moment_2_restrictions_r, moment_2_restrictions_r[,2],  moment_2_restrictions_r[,2]))
    } else {
      moment_3_restrictions_r <- NULL
      moment_4_restrictions_r <- NULL
    }
    
    optim_out <- optim(par = init_val, fn = sem_el_euclid_fit_obj, gr = NULL,
                       y_r = Y, b_r = B, moment_2_restrictions_r = moment_2_restrictions_r,
                       moment_3_restrictions_r = moment_3_restrictions_r, moment_4_restrictions_r = moment_4_restrictions_r,
                       meanEst = meanEst, method = optim.method, control = list(fnscale = -1, reltol = outerTol))
    
    fitted_mod <- sem_el_euclid_fit_weights(optim_out$par, y_r = Y, b_r = B, moment_2_restrictions_r = moment_2_restrictions_r,
                                            moment_3_restrictions_r = moment_3_restrictions_r,
                                            moment_4_restrictions_r = moment_4_restrictions_r, meanEst = meanEst)
    
    p <- c(fitted_mod$p_star)
    
    B.hat <- matrix(0, nrow = V, ncol = V)
    
    if(meanEst){
      mu.hat <- optim_out$par[c(1:V)]
      B.hat[B==1] <- optim_out$par[-c(1:V)]
      temp <- (diag(rep(1,V)) - B.hat) %*% (Y - solve(diag(rep(1,V)) - B.hat, mu.hat))
      
    } else {
      B.hat[B==1] <- optim_out$par
      temp <- (diag(rep(1,V)) - B.hat) %*% Y
    }
    
    
    Omega.hat <-  round(temp %*%diag(p) %*%  t(temp),10)
    lrt <- 1/2 * sum((p * n - 1)^2)
    
    temp <- solve(diag(rep(1,V)) - B.hat)
    
    return(list(B.hat = B.hat,
                Omega.hat = Omega.hat,
                mu.hat = mu.hat,
                p = p,
                converged = !optim_out$convergence & (abs(sum(p)-1) < 1e-8),
                lrt = lrt,
                Sigma.hat = temp %*% Omega.hat %*% t(temp)))
  }
  
  
}