empLike_sandwich <- function(Y, B, Omega, B.hat, p)
{
  beta.index <- getBetaIndex(B)
  omega.index <- getOmegaIndex(Omega)
  V <- dim(B)[1]
  errors <- ((diag(rep(1,V)) - B.hat) %*% Y)  #- mu.hat
  g.i <- rep(0, V + dim(omega.index)[1])

  
  # holds overall diff
  dg <- matrix(0, nrow = V + dim(omega.index)[1],
               ncol =  dim(beta.index)[1])
  g.var <- matrix(0, nrow = V + dim(omega.index)[1],
                  ncol = V + dim(omega.index)[1])
  
  for(n in 1:length(p)){
    g.i[1:V] <- errors[, n]
    g.i[-c(1:V)] <- errors[omega.index[, 1],n] * errors[omega.index[, 2],n]
    
    g.var <- g.var + p[n] * g.i %*% t(g.i)
    
    
    
    dg.i <- matrix(0, nrow = V + dim(omega.index)[1],
                   ncol = dim(beta.index)[1])
    
    # derivative of mean constraints
    for(i in 1:V){
      dg.i[i, which(beta.index[,1] == i)] <- -Y[beta.index[which(beta.index[,1] == i), 2], n]
    }
    
    for(j in 1:dim(omega.index)[1]){
      
      # derivative of beta's
      dg.i[j+V, which(beta.index[,1] == omega.index[j,1])] <- -Y[beta.index[which(beta.index[,1] == omega.index[j,1]), 2], n] *
        errors[omega.index[j,2], n]
      dg.i[j+V, which(beta.index[,1] == omega.index[j,2])] <- -Y[beta.index[which(beta.index[,1] == omega.index[j,2]), 2], n] *
        errors[omega.index[j,1], n]
    }
    dg <- dg + dg.i * p[n]
  }
  
  return(list(dg = dg,
              g.var = g.var,
              sw.var = solve(t(dg) %*% solve(g.var) %*% dg)))
}


#################
getBetaIndex <- function(B)
{
  ret <- NULL
  for(i in 1:dim(B)[1]){
      ret <- rbind(ret, cbind(which(B[,i]==1) , rep(i, sum(B[,i]))))
  }
  return(ret)
}

getOmegaIndex <- function(Omega){ 
  ret <- NULL
  V <- dim(Omega)[1]
  for(i in 1:(V-1)){
    ret <- rbind(ret, cbind(which(Omega[(i+1):V, i]==0)+i, rep(i, sum(Omega[(i+1):V, i]==0))))
  } 
  return(ret)
}

ricf_info <- function(S, BHat, OmegaHat, B, Omega)
{
  V <- dim(B)[1]
  ### Create P Matrix ####
  P <- matrix(0, nrow = V^2, ncol = sum(B))
  temp <- matrix(c(which(c(B)==1),1:sum(B)), nrow = sum(B), ncol = 2)
  P[temp] <- 1
  Q <- matrix(0, nrow = V^2, ncol = (sum(Omega)-V)/2 + V)
  Omega.temp <- matrix(0, nrow = V, ncol = V)
  count <- 0
  for(i in 1:V){
    for(j in i:V){
      if(Omega[i,j] == 1){
        count <- count + 1
        Omega.temp[i,j] <- count
      }
    }
  }
  Omega.temp <- Omega.temp + t(Omega.temp) - diag(diag(Omega.temp))
  Q[cbind(c(1:V^2), c(Omega.temp))] <- 1
  
  
  omega.inv <- solve(OmegaHat)
  eye.minus.b.inv <- solve(diag(rep(1, V)) - BHat)
  Tr <- matrix(0, nrow = V^2, ncol = V ^2)
  Tr[matrix(c(1:V^2, (rep(c(1:V), V)-1) *V + rep(c(1:V), each = V)), ncol = 2)] <- 1
  
  I <- matrix(0, nrow = sum(B) + (sum(Omega)-V)/2 + V,
              ncol = sum(B) + (sum(Omega)-V)/2 + V)
  
  I[1:sum(B), 1:sum(B)] <- t(P) %*% (S %x% omega.inv + (eye.minus.b.inv %x% t(eye.minus.b.inv)) %*% Tr) %*% P
  I[-c(1:sum(B)),-c(1:sum(B))] <- 1/2 * t(Q) %*% (omega.inv %x% omega.inv) %*% Q
  I[c(1:sum(B)),-c(1:sum(B))] <- t(P) %*% ( eye.minus.b.inv %x% omega.inv) %*% Q
  I[-c(1:sum(B)),c(1:sum(B))] <- t(Q) %*% (t(eye.minus.b.inv) %x% omega.inv) %*% P
  
  
  return(I)
}
