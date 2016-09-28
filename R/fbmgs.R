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
fbmgs <- function(Y, moments = 4, meanEst = F, outerTol = 1e-3, innerTol = 1e-3,
                  printStatus = F, fwdNullTol = 30, fwdDirectTol = 15,
                  bwdNullTol = 15, bwdDirectTol = 15){
  V <- dim(Y)[1]
  B <- matrix(0, nrow = V, ncol = V)
  Omega <- diag(rep(1,V))
  
  # hold lr stat of "next best option" from whatever was chosen
  lr.stat <- matrix(0, nrow = V * (V - 1) / 2, ncol = 4)
  
  
  
  ## B and Omega for Pairwise comparisons
  Omega_sub <- diag(rep(1,2))
  B_null <- matrix(c(0, 0,
                     0, 0), nrow = 2, ncol = 2, byrow = T)
  B_i_to_j <- matrix(c(0, 0,
                 1, 0), nrow = 2, ncol = 2, byrow = T)
  B_j_to_i <- t(B_i_to_j)
  
  covarRestrict <- which(Omega_sub == 0 & lower.tri(Omega_sub), arr.ind = T) - 1
  if(moments == 3){
    covarRestrict <- rbind(cbind(covarRestrict, rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,1] ),
                           cbind(covarRestrict, covarRestrict[,2] ))
  } else if (moments == 4) {
    covarRestrict <- rbind(cbind(covarRestrict, rep(-1, dim(covarRestrict)[1]), rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,1] ,rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,2] ,rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,1], covarRestrict[,1]),
                           cbind(covarRestrict, covarRestrict[, 2],covarRestrict[, 2] ))
  }
  
  #### Forward Search ####
  k <- 1
  for(i in 2:V){
    for(j in 1:(i-1)){
      
      if(printStatus){
        cat(paste("\nForward:", i, j))
      }
      
      ## Null Test
      Y.test <- Y[c(i,j), ]
      # edge i <- j
      el_fit_null <- sempl_cd(Y = Y.test, B = B_null, Omega = Omega_sub, meanEst = meanEst,
                           covarRestrict = covarRestrict, tol = innerTol,
                           outerTol = 1e-5)
      
      ## If null model does not fit well, try directed edge
      if(el_fit_null > fwdNullTol){
        
        # edge i <- j
        el_fit_1 <- sempl_cd(Y = Y.test, B = B_i_to_j, Omega = Omega_sub, meanEst = meanEst,
                            covarRestrict = covarRestrict, tol = innerTol,
                            outerTol = 1e-5)
        
        # edge i -> j
        el_fit_2 <- sempl_cd(Y = Y.test, B = B_j_to_i, Omega = Omega_sub, meanEst = meanEst,
                            covarRestrict = covarRestrict, tol = innerTol,
                            outerTol = 1e-5)
        
       if (min(el_fit_1, el_fit_2) <  fwdDirectTol) {
            
            ## Check which edge direction is better
            if(el_fit_1 < el_fit_2){
              B[j, i] <- 1
              lr.stat[k, ] <- c(j, i, el_fit_null - el_fit_1, 0)
              if(printStatus){
                cat(paste("\nAdded directed edge ", i, " -> ", j, "; Best oriented edge: ", round(el_fit_1,3), sep = "" ))
              }
            } else {
              B[i, j] <- 1
              lr.stat[k, ] <- c(i, j, el_fit_null - el_fit_2, 0)
              if(printStatus){
                cat(paste("\nAdded directed edge ", j, " -> ", i, "; Best oriented edge: ", round(el_fit_2,3), sep = "" ))
              }
            }
          } else { 
            ## If we have some dependency that can't be modeled 
            ## by only a directed edge add a bidirected edge
            Omega[i, j] <- Omega[j, i] <- 1
            lr.stat[k, ] <- c(i, j, min(el_fit_1, el_fit_2), 1)
            if(printStatus){
              cat(paste("\nAdded bidirected edge ", j, " <-> ", i, "; Best oriented edge: ", round(min(el_fit_1, el_fit_2), 3), sep = "" )) 
            }
          }
      } else {
        if(printStatus){
          cat(paste("\nNo added edge; Null Model: ", round(el_fit_null,3), sep = "" ))
        }
      }
      
     k <- k + 1
    }
  }
  
  
  lr.stat <- lr.stat[lr.stat[,3] != 0, ]
  lr.stat <- lr.stat[order(lr.stat[,3]),]
  bidirected_prune <- lr.stat[which(lr.stat[,4] == 1), ]
  directed_prune <- lr.stat[which(lr.stat[,4] == 0), ]
  
  
  if(printStatus){
    cat("\nCurrent Graph:\n")
    print(lr.stat)
    cat("\n")
  }
  
  
  
  cat("\n Bidirected Pruning")
  el_fit_saturated <- sempl_cd(Y = Y, B = B, Omega = Omega, meanEst = meanEst, 
                          covarRestrict = .getCovarRestrict(Omega, moments = moments), tol = innerTol,
                          outerTol = outerTol)
  
  ## Knock down bidirected edges to directed edges
  for (z in 1:dim(bidirected_prune)[1]) {
    i <- bidirected_prune[z, 1] 
    j <- bidirected_prune[z, 2] 
    
    
    if (printStatus){
      cat(paste("\nNode ", i ," ", j," : ", sep = ""))
      
    }
    
    Omega.test <- Omega 
    Omega.test[i, j] <- Omega.test[j, i] <- 0
      
    ## Else check if we can remove bidirected edge and add directed edge
    B_i_to_j <- B_j_to_i <- B
    B_i_to_j[j, i] <- 1
    B_j_to_i[i, j] <- 1
        
    el_fit_1 <- sempl_cd(Y = Y, B = B_i_to_j, Omega = Omega.test, meanEst = meanEst, 
                           covarRestrict = .getCovarRestrict(Omega.test, moments = moments), tol = innerTol,
                           outerTol = outerTol)
    
    el_fit_2 <- sempl_cd(Y = Y, B = B_j_to_i, Omega = Omega.test, meanEst = meanEst, 
                           covarRestrict = .getCovarRestrict(Omega.test, moments = moments), tol = innerTol,
                           outerTol = outerTol)
    
    direct.diff <- min(el_fit_1, el_fit_2) - el_fit_saturated
        
    if((min(el_fit_1, el_fit_2) - el_fit_saturated) < bwdDirectTol) {
        ## If completely replacing bidirected edge isn't that bad, then do it 
        Omega <- Omega.test
        if(el_fit_1 < el_fit_2){
          B <- B_i_to_j
          directed_prune <- rbind(directed_prune, c(j, i, min(el_fit_1, el_fit_2) - el_fit_saturated, 0))
          if(printStatus){
            cat(paste("Bidirected changed to directed; ", i, "->" , j, sep = ""))
          }
        } else {
          B <- B_j_to_i
          directed_prune <- rbind(directed_prune, c(i, j, min(el_fit_1, el_fit_2) - el_fit_saturated, 0))
          
          if(printStatus){
            cat(paste("Bidirected changed to directed; ", j, "->" , i, sep = ""))
          }
        }
        if(printStatus){
          cat(paste(" Diff: ", min(el_fit_1, el_fit_2) - el_fit_saturated, sep = ""))
        }
        el_fit_saturated <- min(el_fit_1, el_fit_2)
        
    } else {
      if(printStatus){
        cat(paste("Bidirected edge preserved; ", i, "<->" , j, sep = ""))
        cat(paste(" Diff: ", min(el_fit_1, el_fit_2) - el_fit_saturated, sep = ""))
      }
    }
    
        
    }
    
  

  directed_prune <- directed_prune[order(directed_prune[,3]), ]
  cat("\n Directed Pruning")
  
  for(z in 1:dim(directed_prune)[1]){
    i <- directed_prune[z, 1]
    j <- directed_prune[z, 2]
    
    if (printStatus){
      cat(paste("\nNode ", i ," ", j," : ", sep = ""))
      
    }
    
    B.test <- B
    B.test[i, j] <- 0
    el_fit_null <- sempl_cd(Y = Y, B = B.test, Omega = Omega, meanEst = meanEst, 
                        covarRestrict = .getCovarRestrict(Omega, moments = moments), tol = innerTol,
                        outerTol = outerTol)
      
    null.diff <- el_fit_null - el_fit_saturated

      
    if(null.diff < bwdNullTol) {
        el_fit_saturated <- el_fit_null
        B <- B.test
        if(printStatus){
          cat(paste("Directed edge removed; ", i, " ", j, sep = ""))
          cat(paste(" Diff: ", null.diff, sep = ""))
        }
    } else {
      if(printStatus){
        cat(paste("Directed edge preserved; ", j, "->", i, sep = ""))
        cat(paste(" Diff: ", null.diff, sep = ""))
      }
    }
  }
  
  return(list(B = B, Omega = Omega))
}









.getCovarRestrict <- function(Omega, moments = 4){
  covarRestrict <- which(Omega == 0 & lower.tri(Omega), arr.ind = T) - 1
  if(moments == 3){
    covarRestrict <- rbind(cbind(covarRestrict, rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,1, drop = F] ),
                           cbind(covarRestrict, covarRestrict[,2, drop = F] ))
  } else if (moments == 4) {
    covarRestrict <- rbind(cbind(covarRestrict, rep(-1, dim(covarRestrict)[1]), rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,1, drop = F] ,rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,2, drop = F] ,rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,c(1,1), drop = F]),
                           cbind(covarRestrict, covarRestrict[,c(2,2), drop = F] ))
  }
  return(covarRestrict)
}