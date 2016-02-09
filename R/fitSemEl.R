#' Fitting Linear Structural Equation Models
#' 
#' 
#' Fits Empirical Likelihood estimates for given SEM
#' 
#'  
#' @param B V by V matrix with {0,1} giving structure of directed edges
#' @param Omega V by V matrix with {0,1} giving structure of bi-directed edges
#' @param B.hat V by V matrix giving initial edges weights for directed edges. If BInit is NULL,
#'    a default initialization will be used. 
#' @param Omega.hat V by V matrix giving initial edge weights for bi-directed edges
#' @param mu.hat Initial value for means  
#' @param optim.method
#' @param tol
#' @param maxInnerIter
#' @param maxOuterIter
#' @param meanEst
#' @param type profile, naive, euclid
#' @param high_moments
#' @export
fitEL <- function(Y, B, Omega, B.hat = NULL, Omega.hat = NULL, mu.hat = NULL, optim.method = "BFGS",
                  tol = 1e-6, maxInnerIter = 500, maxOuterIter = 100, meanEst = 0, type = "profile" , high_moments = 2){
  
  
  
  if(type == "profile"){
    V <- dim(B)[1]
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
    
    
    optim_out <- optim(par = init_val, fn = sem_el_fit_obj, gr = NULL,
                       y_r = Y, b_r = B, moment_2_restrictions_r = moment_2_restrictions_r,
                       moment_3_restrictions_r = moment_3_restrictions_r, moment_4_restrictions_r = moment_4_restrictions_r,
                       tol = tol, max_iter = maxInnerIter, meanEst = meanEst, method = optim.method, control = list(fnscale = -1, trace = 0))
    
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
      B.hat[B==1] <- optim_out$par
      temp <- (diag(rep(1,V)) - B.hat) %*% Y  %*% diag(sqrt(p))
    }
    
    
    
    
    Omega.hat <-  round(temp %*% t(temp),10)
    lrt <- -2 * sum(log(p * n))
    
    temp <- solve(diag(rep(1,V)) - B.hat)
    
    return(list(B.hat = B.hat,
                Omega.hat = Omega.hat,
                mu.hat = mu.hat,
                p = p,
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
                       y_r = Y, b_r = B, moment_2_restrictions_r = moment_2_restrictions_r,
                       moment_3_restrictions_r = moment_3_restrictions_r, moment_4_restrictions_r = moment_4_restrictions_r,
                       tol = tol, max_Iter = maxInnerIter, meanEst = meanEst, method = optim.method, control = list(fnscale = -1, trace = 0))
    
#     fitted_mod <- sem_el_naive_fit_weights(optim_out$par, y_r = Y, b_r = B, moment_2_restrictions_r = moment_2_restrictions_r,
#                                      moment_3_restrictions_r = moment_3_restrictions_r, moment_4_restrictions_r = moment_4_restrictions_r,
#                                      meanEst = meanEst, tol = tol, max_Iter = maxInnerIter)
    
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
      out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                            OmegaInit = NULL, sigConv = 0, maxIter = 0,
                            msgs = FALSE, omegaInitScale = .9)
      init_val = c(out_ricf_init$BHat[B==1])
    } else {
      init_val = c(B.hat[B==1])
    }
    
    # For the mean (if needed)
    if(meanEst & !is.null(mu.hat)){
      init_val <- c(mu.hat, init_val)
    } else if(meanEst){
      init_val <- c(rowMeans(Y), init_val)
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
                       meanEst = meanEst, method = optim.method, control = list(fnscale = -1, trace = 0))
    
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





# # Calculate weights at a specified point
# fitEL.weights <- function(Y, B, Omega, B.hat, mu.hat = NULL,
#                           tol = 1e-6, maxInnerIter = 100){
#   
#   num_dual_vars = V + sum(Omega == 0)/2
#   
#   val <- B.hat[which(B==1)]
#   
#   fitted_mod <- sem_el_fit_weights(val, y_r = Y, omega_r = Omega, b_r = B,
#                                    dual_r = rep(0, num_dual_vars), tol = tol, max_iter = maxInnerIter, meanEst = !is.null(mu.hat))
#   
#   p <- c(1/fitted_mod$d)
#   temp <- (diag(rep(1,V)) - B.hat) %*% Y %*% diag(sqrt(p))
#   Omega.hat <-  round(temp %*% t(temp),10)
#   lrt <- -2 * sum(log(p * n))
#   
#   return(list(Omega.hat = Omega.hat,
#               p = p,
#               converged = (abs(sum(p)-1) < 1e-8),
#               lrt = lrt))
#   
# }