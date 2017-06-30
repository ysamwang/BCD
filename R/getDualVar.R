getDualVar <- function(Y, B, Omega, mod)
{
  # Index of non-zero B's in proper order (column by column)
  beta.index <- .getBetaIndex(B)
  # Index of non-zero Omega's
  omega.index <- .getOmegaIndex(Omega)
  # number of Nodes
  V <- dim(B)[1]
  # errors for each node and observation
  errors <- ((diag(rep(1,V)) - mod$B.hat) %*% Y)  #- mu.hat
  p <- mod$p
  
  # pre-allocate memory (note that dg.i is cleared out in each iteration of the loop)
  dg <- matrix(0, nrow = V + V * (V + 1) /2,
               ncol =  dim(beta.index)[1] + dim(omega.index)[1])
  g.var <- matrix(0, nrow = V + V * (V + 1) /2,
                  ncol = V + V * (V + 1) /2)
  g.i <- rep(0, V + V * (V+1)/2)
  
  
  
  # indexes of non-replicated entries of Omega
  # This includes 0 entries as well whereas Omega.index only includes non-zero entries
  covar.temp <- matrix(0, ncol = 2, nrow = V * (V + 1) /2)
  count <- 1
  for(i in 1:V){
    # TO-DO: modify code so it is similar to Omega.index in cbinding to null
    # probably more readable (could be tad slower)
    covar.temp[c(count:(count + V -i)), c(1:2)] <- cbind(c(i:V), rep(i, (V-i+1)))
    count <- count + (V-i+1)
  }
  
  
  for(n in 1:length(p)){
    ### Form m(X_i, \theta ) ### 
    
    # mean restrictions
    g.i[1:V] <- errors[, n]
    # covariance restrictions
    g.i[-c(1:V)] <- errors[covar.temp[, 1], n] * errors[covar.temp[, 2] ,n] - mod$Omega.hat[covar.temp]
    
    g.var <- g.var + p[n] * g.i %*% t(g.i)
    
    
    ### Now form dm/dtheta
    dg.i <- matrix(0, nrow = V + V * (V + 1)/2,
                   ncol = dim(beta.index)[1] + dim(omega.index)[1])
    
    # derivative of mean constraints
    for(k in 1:V){
      dg.i[k, which(beta.index[, 1] == k)] <- -Y[beta.index[which(beta.index[, 1] == k), 2], n]
    }
    
    # derivative of covariance constraints
    for(j in 1:(V * (V+1) / 2)){
      dg.i[V + j, which(beta.index[,1] == covar.temp[j,1])] <- -Y[beta.index[which(beta.index[,1] == covar.temp[j, 1]), 2], n] * errors[covar.temp[j,2], n]
      dg.i[V + j, which(beta.index[,1] == covar.temp[j,2])] <- -Y[beta.index[which(beta.index[,1] == covar.temp[j, 2]), 2], n] * errors[covar.temp[j,1], n]
      if(Omega[covar.temp[j,1],covar.temp[j,2]] == 1){
        dg.i[V + j, dim(beta.index)[1] + which(covar.temp[j,1] == omega.index[,1] & covar.temp[j,2] == omega.index[,2])] <- -1 
      }
    }
    dg <- dg + dg.i * p[n]
  }
  V1 <- solve(t(dg) %*% solve(g.var) %*% dg)
  U <- solve(g.var) %*% (diag(rep(1, dim(g.var)[1]))  -  dg %*% V1 %*% t(dg) %*% solve(g.var))
  return(list(U = U, V = V1, omegaIndex = omega.index, covarTemp = covar.temp))
}