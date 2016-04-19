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
forwardBackwardFull <- function(Y, moments = 3, meanEst = 1,
                               optimTol = 1e-4, innerTol = 1e-4,
                               printStatus = T, exch = F, vsNull = T,
                               fwdTol = 30, fwdType = "profile", bwdTol = 15, bwdType = "euclid"){
  
  # number of variables
  V <- dim(Y)[1]
  # Initial estimates of B and Omega
  B <- matrix(0, nrow = V, ncol = V)
  Omega <- matrix(1, nrow = V, ncol = V)
  
  # list of node pairs 
  lr.stat <- cbind(t(combn(1:V, 2)), rep(0, choose(V, 2)))
  
  stoppingCrit <- fwdTol + 1
  while(dim(lr.stat)[1] > 0 & stoppingCrit > fwdTol) {
    cat("\n=== Forward Step (")
    cat(sum(B))
    cat(" edges added) ===")
    
    # Run through all available edges and see which one to add
    for(k in 1:dim(lr.stat)[1]){
      i <- lr.stat[k, 1]
      j <- lr.stat[k, 2]
      
      if(printStatus){
        cat(paste("\n\tChecking:", i, j))
      }
      
      B1 <- B2 <- B3 <- B
      # add edges for testing
      B1[i, j] <- B2[j, i] <- 1
      # punch out holes in omega
      Omega[i, j] <- Omega[j, i] <- 0
      ### Tests ####
      
      # edge i <- j
      el_fit1 <- fitEL(Y = Y, B = B1, Omega = Omega,
                       meanEst = meanEst, type = fwdType,
                       high_moments = moments, optim.method = "BFGS",
                       tol = innerTol, outerTol = optimTol)
      
      # edge i -> j
      el_fit2 <- fitEL(Y = Y, B = B2, Omega = Omega,
                       meanEst = meanEst, type = fwdType,
                       high_moments = moments, optim.method = "BFGS",
                       tol = innerTol, outerTol = optimTol)
      
      
        
        if (printStatus){
          cat(paste("\n\tWith ",i, "," ,j, ": ", round(el_fit1$lrt,3), sep = ""))
          cat(paste("\tWith ",j, "," ,i, ": ", round(el_fit2$lrt,3), sep = ""))
        }
        
        
        if (el_fit2$lrt >  el_fit1$lrt  ){
          # B[i, j] <- 1
          lr.stat[k, ] <- c(i, j, el_fit2$lrt -  el_fit1$lrt)
        } else {
          # B[j, i] <- 1
          lr.stat[k, ] <- c(j, i, el_fit1$lrt -  el_fit2$lrt)
        }
      
      Omega[i, j] <- Omega[j,i] <- 1
    }
    
    # If there is one that has a very large discrepancy, add it
    stoppingCrit <- max(lr.stat[,3])
    if(stoppingCrit > fwdTol){
      update <- which.max(lr.stat[,3])
      B[lr.stat[update, 1], lr.stat[update, 2] ] <- 1
      
      cat("\nAdded: ")
      cat(paste(lr.stat[update, 1], ", ",lr.stat[update, 2] ,sep = "")) 
      cat("\nDiff: ")
      print(lr.stat[update, 3])
      

      
      Omega[lr.stat[update, 1], lr.stat[update, 2]] <- 
      Omega[lr.stat[update, 2], lr.stat[update, 1]] <- 0
      
      lr.stat <- lr.stat[-update, , drop = FALSE]
    }
  }
  
  # Done with Fwd Growing
  cat("\n=== Done with Fwd ===\n")
  
  
  # Start Backward Pruning
  print(B)
  cat("\n=== Begining Bwd ===\n")
  
  
  #### Backward Prune ####


  Omega <- diag(rep(1, V))
  stoppingCrit <- bwdTol - 1
  
  # Current Edges
  lr.stat <- cbind(which(B == 1, arr.ind = T), rep(0, sum(B)))
  
  # Fit Model
  el_fit2 <- fitEL(Y = Y, B = B, Omega = Omega,
                   meanEst = meanEst, type = fwdType,
                   high_moments = moments, optim.method = "BFGS",
                   tol = innerTol, outerTol = optimTol)
  full.lrt <- el_fit2$lrt

  while(dim(lr.stat)[1] > 0 & stoppingCrit < bwdTol) {
    cat("\n=== Backward Step (")
    cat(sum(B))
    cat(" edges remaining) ===")
    cat("\nCurrent LR: ")
    cat(round(full.lrt, 3))
    cat("\n")
    
    
    # Test each edge
    for(k in 1:dim(lr.stat)[1]){
      
      i <- lr.stat[k, 1]
      j <- lr.stat[k, 2]
      
      if(printStatus){
        cat(paste("\n\tChecking:", i, j))
      }
      
      
      # See if we need to fit this model or not
      if( lr.stat[k, 3] - full.lrt < bwdTol ){
        B1 <- B
        # Remove edge for testing
        B1[i, j] <- 0
          
        el_fit1 <- fitEL(Y = Y, B = B1, Omega = Omega,
                         meanEst = meanEst, type = fwdType,
                         high_moments = moments, optim.method = "BFGS",
                         tol = innerTol, outerTol = optimTol)
        
        
        if (printStatus) {
          cat(paste("\n\twithout: ", round(el_fit1$lrt,3), sep = ""))
          
          }
        # Record lr stat
        lr.stat[k, ] <- c(i, j, el_fit1$lrt)
        
        } else { # Since model is nested in previous model, if it's already too significant, skip
        if (printStatus){
          cat(paste("\n\twithout (no-refit): ", round(lr.stat[k, 3],3), sep = ""))
        }
      }
    } # Done testing each edge
    
    # check if smallest lr-stat is within full model by certain tolerance
    stoppingCrit <- min(lr.stat[,3]) - full.lrt
    
    if(stoppingCrit < bwdTol){
      
      # get index of smallest lr-stat
      update <- which.min(lr.stat[,3])
      # Update B
      B[lr.stat[update, 1], lr.stat[update, 2] ] <- 0
      
      if(printStatus){
        cat("\nEdge removed: ")
        print(lr.stat[update, 1:2])
        cat("\nDiff: ")
        print(lr.stat[update, 3] - full.lrt)
      }
      
      # record lr stat of model (so we don't have to refit)
      full.lrt <- lr.stat[update, 3]
      # delete edge from list of future edges to check
      lr.stat <- lr.stat[-update, , drop = FALSE]      
    }
  }
  
  return(B)
}
  
  
