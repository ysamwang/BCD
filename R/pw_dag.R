sempl_pwdag_search <- function(Y, moments = 3,
                            null_cutoff = qchisq(.9, df = 2 + 1 + moments * 2),
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
        
        if(out.12$lr < out.21$lr ) {
          
          ## Set edge in main B
          B[E_remaining[i, 1], E_remaining[i, 2]] <- 1
          
          ## note that a change has been made
          change_made = T
          null_edges <- c(null_edges, i)
        } else  {
          
          ## Set edge in main B
          if(printStatus){
            print(paste("Setting edge:", E_remaining[i , 1], "->", E_remaining[i , 2], collapse = " "))
          }
          B[E_remaining[i, 2], E_remaining[i , 1]] <- 1
          
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
  
  return( B)
}




