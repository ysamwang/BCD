forwardBackwardDiscoveryEuclid <- function(Y, moments = 3, tol = 30, type = "euclid", sat = T){
  V <- dim(Y)[1]
  B <- matrix(0, nrow = V, ncol = V)
  if(sat){
    Omega <- matrix(1, nrow = V, ncol = V)
  } else { 
    Omega <- diag(rep(1,V))
  }
  
  # forward selection
  for(i in 1:(V-1)){
    for(j in (i+1):V){
      cat(paste("\nForward:", i,j))
      B1 <- B
      B1[i, j] <- 1
      
      if(sat){
        Omega[i,j] <- Omega[j,i] <- 0
      }
      
      el_fit1 <- fitEL(Y = Y, B = B1, Omega = Omega,
                       meanEst = 1, type = type, high_moments = moments, optim.method = "BFGS",
                       tol = 1e-4, outerTol = 1e-2)
      
      B2 <- B
      B2[j,i] <- 1
      el_fit2 <- fitEL(Y = Y, B = B2, Omega = Omega,
                       meanEst = 1, type = type, high_moments = moments, optim.method = "BFGS", tol = 1e-4,
                       outerTol = 1e-2)
      
      el_fit3 <- fitEL(Y = Y, B = B, Omega = Omega,
                       meanEst = 1, type = type, high_moments = moments, optim.method = "BFGS", tol = 1e-4,
                       outerTol = 1e-2)
      
      
      diff <- el_fit2$lrt - el_fit1$lrt
      cat(paste("\n\t",i, "," ,j, ": ", round(el_fit1$lrt,3), sep = ""))
      cat(paste("\n\t",j, "," ,i, ": ", round(el_fit2$lrt,3), sep = ""))
      cat(paste("\n\t",i, "," ,j, " (none) : ", round(el_fit3$lrt,3), sep = ""))
      if( (el_fit3$lrt - min(el_fit2$lrt - el_fit1$lrt)) > tol ){
        if(el_fit2$lrt > el_fit1$lrt){
          B <- B1
        } else {
          B <- B2
        }
      }
      
      if(sat){
        Omega[i,j] <- Omega[j,i] <- 1
      }
    }
  }
  
  cat("\n")
  # Backward Pruning
  # forward selection
  
  Omega <- diag(rep(1,V))
  for(i in 1:(V-1)){
    for(j in (i+1):V){
      
      if(B[i,j] + B[j,i] == 1){
        cat(paste("Backward:", i,j))
        B1 <- B
        B1[i, j] <- B1[j, i] <- 0
        el_fit1 <- fitEL(Y = Y, B = B1, Omega = Omega,
                         meanEst = 1, type = type, high_moments = moments, optim.method = "BFGS", tol = 1e-4,
                         outerTol = 1e-4)
        
        B2 <- B
        el_fit2 <- fitEL(Y = Y, B = B2, Omega = Omega,
                         meanEst = 1, type = type, high_moments = moments, optim.method = "BFGS", tol = 1e-4, outerTol = 1e-4)
        
        
        diff <- el_fit1$lrt - el_fit2$lrt
        cat(paste(" ", round(diff,3), "\n"))
        if( diff < tol){
          B <- B1
        }
        
      }
    }
  }
  
  return(B)
  
  
}