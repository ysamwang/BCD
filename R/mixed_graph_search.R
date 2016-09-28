sempl_mg_search <- function(Y, moments = 3,
                            null_cutoff = qchisq(.95, df = 2 + 1 + (moments-1) * 2),
                            t_cutoff = 3, printStatus = F){
  V <- dim(Y)[1]
  B <- matrix(0, nrow = V, ncol = V)
  
  
  
  ## Scan through until all directed edges are identified
  E_remaining <- which(lower.tri(B), arr.ind = T)
  change_made <- T
  B.12 <- B.21 <- B.00 <-  matrix(0, nrow = 2, ncol = 2)
  B.12[1, 2] <- 1
  B.21[2, 1] <- 1
  Omega.diag <- diag(rep(1,2))
  covarRestrict <- getCovarRestrict(Omega.diag, moments = moments)
  
  Y.current <- Y
  
  while(change_made & (dim(E_remaining)[1] > 0)){

    change_made = F
    null_edges <- c()
    
    ## Go through each edge pairing
    for(i in 1:dim(E_remaining)[1]){
      
      ## Fit competing directed models
      Y.test <- Y.current[E_remaining[i, ], ]
      
      out.00 <- sempl(Y.test, B.00, Omega.diag, meanEst = F,
                         covarRestrict = covarRestrict) 
      
      if(printStatus){
        print(paste(E_remaining[i, 1], " -/-", E_remaining[i, 2], " : ", round(out.00$lr, 3)), collapse = "")  
      }
      
      if( out.00$lr < null_cutoff ) {
        ## Null model passes
        null_edges <- c(null_edges, i)
      } else {
        
        out.12 <- sempl(Y.test, B.12, Omega.diag, meanEst = F,
                        covarRestrict = covarRestrict)
        
        out.21 <- sempl(Y.test, B.21, Omega.diag, meanEst = F,
                        covarRestrict = covarRestrict)
        if(printStatus){
          print(paste(E_remaining[i, 1], " ->", E_remaining[i, 2], " : ", round(out.21$lr, 3)), collapse = "")  
          print(paste(E_remaining[i, 2], " ->", E_remaining[i, 1], " : ", round(out.12$lr, 3)), collapse = "")
        }
        
        if(out.12$lr < null_cutoff & out.21$lr > null_cutoff) {
            
            ## Set edge in main B
            B[E_remaining[i, 1], E_remaining[i, 2]] <- 1
            
            ## Update error estimates
            Y.current[c(E_remaining[i, 1], E_remaining[i, 2]), ] <- Y.test - out.12$B.hat %*% Y.test
            
            ## note that a change has been made
            change_made = T
            null_edges <- c(null_edges, i)
            if(printStatus){
              print(paste("Setting edge:", E_remaining[i , 2], "->",E_remaining[i , 1], collapse = " "))
            }
          } else if(out.21$lr < null_cutoff & out.12$lr > null_cutoff) {
            
            ## Set edge in main B
            if(printStatus){
              print(paste("Setting edge:", E_remaining[i , 1], "->", E_remaining[i , 2], collapse = " "))
            }
            B[E_remaining[i, 2], E_remaining[i , 1]] <- 1
            
            ## Update error estimates
            Y.current[c(E_remaining[i, 1], E_remaining[i, 2]), ] <- Y.test - out.21$B.hat %*% Y.test
            
            ## note that an edge has changed
            change_made = T
            null_edges <- c(null_edges, i)
            
          }
        
        }
      
      
    }
    
    
    if(length(null_edges) > 0){
      E_remaining <- E_remaining[-null_edges, , drop = F]
    }
    if(printStatus){
      print(E_remaining)
    }
  }
  
  ## All leftover edges must be bidirected
  Omega <- matrix(0, nrow = V, ncol = V)
  Omega[E_remaining] <- 1
  Omega <- Omega + t(Omega) + diag(rep(1,V))
  

  return(list( B = B, Omega = Omega))
  
  full_model <- sempl(Y, B, Omega, meanEst = F)
  
  v <- var.el(Y, B, Omega, full_model$B.hat, full_model$Omega.hat, p = full_model$p)
  
  t.val <- abs(c(full_model$B.hat[which(B == 1)],
                 full_model$Omega.hat[which(lower.tri(Omega, diag = T) & Omega == 1)] )) /
                 sqrt(diag(v))
  
  if(printStatus){
    print(t.val)
  }
  
  B.ret <- matrix(0, nrow = V, ncol = V) 
  B.ret[ which(B == 1)[which(t.val[1:sum(B)] > t_cutoff)] ] <- 1
  
  if(printStatus){
    cat("Pruning directed edges:\n")
    print(which(B == 1, arr.ind = T)[which(t.val[1:sum(B)] > t_cutoff), ] )
  }
  
  Omega.ret <- matrix(0, nrow = V, ncol = V)
  Omega.ret[which(Omega == 1 & lower.tri(Omega, diag = T) )[which(t.val[-c(1:sum(B))] > t_cutoff)]] <- 1
  
  if(printStatus){
    cat("Pruning bidirected edges:\n")
    print(which(Omega == 1 & lower.tri(Omega, diag = T) , arr.ind = T)[which(t.val[-c(1:sum(B))] > t_cutoff), ])
  }
  
  
  
  ## Set diagonal elements to 1
  Omega.ret <- Omega.ret + t(Omega.ret)
  diag(Omega.ret) <- 1
  return(list( B = B.ret, Omega = Omega.ret))
}

  
  

