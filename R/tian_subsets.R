tian.subsets <- function(Y, B, Omega, Omega.hat){
  V <- dim(Y)[1]

  omega.inv <- solve(Omega.hat)
  remaining.nodes <- c(1:V)
  
  
  # Get connected components and parents of connected components
  connected.subsets <- list()
  while(length(remaining.nodes) > 0 ){
    i <- remaining.nodes[1]
    connected.i <- which(omega.inv[i,] != 0)
    if(length(connected.i) == 1){
      pa.con.i <- which(B[connected.i,]==1)
    } else {
      pa.con.i <- unique(unlist(apply(B[connected.i,], MAR = 1, function(x){which(x==1)})))
    }
    connected.subsets[[length(connected.subsets) + 1]] <- list(set = connected.i, parents = setdiff(unique(pa.con.i), connected.i))
    remaining.nodes <- setdiff(remaining.nodes, connected.i)
  }
  B.est <- Omega.est <- matrix(0, nrow = V, ncol =  V)
  
  for(i in 1:length(connected.subsets))
  {
    set <- connected.subsets[[i]]$set
    pa <- connected.subsets[[i]]$parents
    B.mod <- B[c(set, pa), c(set, pa) ]
    B.mod[-c(1:length(set)), ] <- 0 
    Y.mod <- Y[c(set, pa), ]  
    Omega.mod <- Omega[c(set, pa), c(set, pa) ] 
    Omega.mod[-c(1:length(set)), ] <- 1
    Omega.mod[,-c(1:length(set))] <- 1
    mod <- fitEL(Y.mod, B.mod, Omega.mod)
    B.est[set, c(set,pa)]  <- mod$B.hat[c(1:length(set)),]
    Omega.est[set, set]  <- mod$Omega.hat[c(1:length(set)), c(1:length(set))]
  }
  
  return(list(B.hat = B.est, Omega.hat = Omega.est))
}


