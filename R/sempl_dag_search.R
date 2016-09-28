sempl_dag_search <- function(Y, moments = 3,
                             null_cutoff = qchisq(.9, df = 2 + 1 + moments * 2),
                             t_cutoff = 3, printStatus = F){
  ordering <- c()
  V <- dim(Y)[1]
  B <- matrix(0, nrow = V, ncol = V)
  
  remaining.edges <- c(1:V)
  
  B.12 <- B.21 <- B.00 <-  matrix(0, nrow = 2, ncol = 2)
  B.12[1, 2] <- 1
  B.21[2, 1] <- 1
  Omega.diag <- diag(rep(1,2))
  covarRestrict <- getCovarRestrict(Omega.diag, moments = moments)
  
  Y.current <- Y
  
  while(length(remaining.edges) > 1){
    
    change_made = F
    null_edges <- c()
    
    lrt <- B.temp <-  matrix(0, nrow = length(remaining.edges), ncol = length(remaining.edges))
    colnames(lrt) <- rownames(lrt) <- remaining.edges
    
    ## Go through each edge pairing
    for(i in remaining.edges){
      for(j in remaining.edges){
        if(i != j){
          Y.test <- Y.current[c(i,j), ]  
          
          
          out.12 <- sempl(Y.test, B.12, Omega.diag, meanEst = F,
                          covarRestrict = covarRestrict)
          
          i.index <- which(remaining.edges == i)
          j.index <- which(remaining.edges == j)
          
          lrt[i.index, j.index] <- out.12$lr
          
          B.temp[i.index, j.index] <- out.12$B.hat[1,2]
        }
      }
    }
    
    root <- which.min(apply(lrt, MAR = 2, max))
    if(printStatus){
      print(lrt)
      cat(paste("Root Found: ", remaining.edges[root], sep = ""))
      cat(paste("; LRT: ", round(min(apply(lrt, MAR = 2, max)), 4) , sep = ""))
      cat("\n")
    }
    # cat(round(max(lrt[,root]),3))
    
    B[remaining.edges[-root], remaining.edges[root]] <- 1
    ordering <- c(ordering, remaining.edges[root])
    
    Y.current[remaining.edges, ] <- Y.current[remaining.edges, ,drop = F] - B.temp[, root, drop = F] %*% Y.current[remaining.edges[root], , drop = F]
    remaining.edges <- remaining.edges[-root]
  }
  ordering <- c(ordering, remaining.edges)
  
  
  Omega <- diag(rep(1, V))
  full_model <- sempl(Y, B, Omega, meanEst = F)
  
  v <- NULL
  
  try(v <- var.el(Y, B, Omega, full_model$B.hat, full_model$Omega.hat, p = full_model$p))
  
  if( !is.null(v) ){
    t.val <- abs(full_model$B.hat[which(B == 1)]) / sqrt(diag(v))[1:sum(B)]
    
    B.ret <- matrix(0, nrow = V, ncol = V) 
    B.ret[ which(B == 1)[which(t.val > t_cutoff)] ] <- 1
  } else {
    B.ret <- B
  }
  return(list(B = B.ret, order = ordering))
}




