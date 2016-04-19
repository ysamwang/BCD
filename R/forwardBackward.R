#' Discovery of Linear Structural Equation Models
#' 
#' Produces estimated SEM from edge-wise tests
#'  
#' @param Y data in V x N matrix
#' @param moments number of moments to test in estimating equations
#' @param fwdTol tolerance for difference on forward search
#' @param bwdTol tolerance for difference on backward search
#' @param fwdType type of fit for forward search: profile or euclid
#' @param bwdType type of fit for backward search: profile or euclid
#' @param meanEst whether or not to fit a mean
#' @param optimTol tolerance for optim
#' @param innerTol tolerance for dual problem
#' @param printStatus whether or not to print output
#' @export
forwardBackwardDiscovery <- function(Y, moments = 3, meanEst = 1,
                                     optimTol = 1e-4, innerTol = 1e-4,
                                     printStatus = T, exch = F, vsNull = T,
                                     fwdTol = 30, fwdType = "profile", bwdTol = 15, bwdType = "euclid"){
  V <- dim(Y)[1]
  B <- matrix(0, nrow = V, ncol = V)
  
  # hold calculated differences 
  lr.stat <- matrix(0, nrow = V, ncol = V)
  
  if(exch){
    Omega <- diag(rep(1,2))
    B1 <- matrix(c(0, 1,
                   0, 0), nrow = 2, ncol = 2, byrow = T)
    B2 <- t(B1)
    B3 <- matrix(0, nrow = 2, ncol = 2)
  } else {
    Omega <- matrix(1, nrow = V, ncol = V)
    Y.test <- Y
  }
  
  #### Forward Search ####
  
    for(i in 2:V){
      for(j in 1:(i-1)){
        
        if(printStatus){
          cat(paste("\nForward:", i, j))
        }
        
        
        # set baseline B may include previous edges or not
        if(exch){
          Y.test <- Y[c(i,j), ]
        } else {
          B1 <- B2 <- B3 <- B
          # add edges for testing
          B1[i, j] <- B2[j, i] <- 1
          # punch out holes in omega
          Omega[i, j] <- Omega[j, i] <- 0
        }
        
        
        
        
        
        
        ### Tests ####
        
        # edge i <- j
        el_fit1 <- fitEL(Y = Y.test, B = B1, Omega = Omega,
                         meanEst = meanEst, type = fwdType,
                         high_moments = moments, optim.method = "BFGS",
                         tol = innerTol, outerTol = optimTol)
        
        # edge i -> j
        el_fit2 <- fitEL(Y = Y.test, B = B2, Omega = Omega,
                         meanEst = meanEst, type = fwdType,
                         high_moments = moments, optim.method = "BFGS",
                         tol = innerTol, outerTol = optimTol)
        
        
        if(vsNull){
          # no edge i -\- j
          el_fit3 <- fitEL(Y = Y.test, B = B3, Omega = Omega,
                           meanEst = meanEst, type = fwdType,
                           high_moments = moments, optim.method = "BFGS",
                           tol = innerTol, outerTol = optimTol)
          
          if (printStatus){
            cat(paste("\n\t",i, "," ,j, ": ", round(el_fit1$lrt,3), sep = ""))
            cat(paste("\n\t",j, "," ,i, ": ", round(el_fit2$lrt,3), sep = ""))
            cat(paste("\n\t",i, "," ,j, " (none) : ", round(el_fit3$lrt,3), sep = ""))
          }
          
          
          # update B if adding edge decreases LR stat by some tolerance
          if( (el_fit3$lrt - min(el_fit2$lrt, el_fit1$lrt)) > fwdTol ){
            
            # add either edge i <- j or i -> j 
            if(el_fit2$lrt > el_fit1$lrt){
              B[i, j] <- 1
              lr.stat[i, j] <- (el_fit3$lrt - min(el_fit2$lrt, el_fit1$lrt))
            } else {
              B[j, i] <- 1
              lr.stat[j, i] <- (el_fit3$lrt - min(el_fit2$lrt, el_fit1$lrt))
            }
          }
          
          
        } else {
          
          if (printStatus){
            cat(paste("\n\t",i, "," ,j, ": ", round(el_fit1$lrt,3), sep = ""))
            cat(paste("\n\t",j, "," ,i, ": ", round(el_fit2$lrt,3), sep = ""))
          }
          
          
          if (el_fit2$lrt -  el_fit1$lrt  > fwdTol){
            B[i, j] <- 1
            lr.stat[i, j] <- el_fit2$lrt -  el_fit1$lrt
          } else if (el_fit1$lrt -  el_fit2$lrt  > fwdTol) {
            B[j, i] <- 1
            lr.stat[j, i] <- el_fit1$lrt -  el_fit2$lrt
          }
      }
        
      }
    }
  
  
  
  order_of_prune <- cbind(which(lr.stat != 0 , arr.ind = T), lr.stat[which(lr.stat != 0)])
  order_of_prune <- order_of_prune[order(order_of_prune[,3]),]
  
  if(printStatus){
    cat("\nCurrent B:\n")
    print(order_of_prune)
    cat("\n")
  }
  
  
  Omega <- diag(rep(1, V))
  # Backward Pruning
  
  for(z in 1:length(order_of_prune[,1])){  
        i <- order_of_prune[z, 1]
        j <- order_of_prune[z, 2]
        
        B1 <- B
        B1[i, j] <- 0
        
        el_fit1 <- fitEL(Y = Y, B = B1, Omega = Omega,
                         meanEst = meanEst, type = bwdType,
                         high_moments = moments, optim.method = "BFGS",
                         tol = innerTol, outerTol = optimTol)        
        B2 <- B
        el_fit2 <- fitEL(Y = Y, B = B2, Omega = Omega,
                             meanEst = meanEst, type = bwdType,
                             high_moments = moments, optim.method = "BFGS",
                             tol = innerTol, outerTol = optimTol)
            
          diff <- el_fit1$lrt - el_fit2$lrt
            
            
            if(printStatus){
              cat(paste("Backward:", i, j))
              cat(paste(" ", round(diff,3), "\n"))
            }
            
            if( diff < bwdTol){
              B <- B1
            }
            
          }
          
          return(B)
}