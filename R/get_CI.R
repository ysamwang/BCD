#' Empirical Likelihood Methods for Linear Structural Equation Models
#' 
#' Calculate assymptotic variance for non-structural zeros in B and Omega matrix. The ordering of the estimates
#' in the variance matrix follows standard vec(-) convention indexing over rows first then columns. The elements of B
#' are ordered before the elements of Omega. Duplicate elements of Omega are not included.
#'
#' @param Y V by N matrix with observed data    
#' @param B V by V matrix with {0,1} giving structure of directed edges
#' @param Omega V by V matrix with {0,1} giving structure of bi-directed edges
#' @param B.hat V by V matrix giving estimated edges weights for directed edges
#' @param Omega.hat V by V matrix giving estimated edge weights for bi-directed edges
#' @param p vector of length N containing the estimated empirical likelihood weights
#' @return The asymptotic covariance matrix as derived in Qin and Lawless 1994
#' @export
getELVariance <- function(Y, B, Omega, B.hat, Omega.hat, p)
{
  # Index of non-zero B's in proper order (column by column)
  beta.index <- getBetaIndex(B)
  # Index of non-zero Omega's
  omega.index <- getOmegaIndex(Omega)
  # number of Nodes
  V <- dim(B)[1]
  # errors for each node and observation
  errors <- ((diag(rep(1,V)) - B.hat) %*% Y)  #- mu.hat
  
  
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
    g.i[-c(1:V)] <- errors[covar.temp[, 1], n] * errors[covar.temp[, 2] ,n] - Omega.hat[covar.temp]
    
    g.var <- g.var + p[n] * g.i %*% t(g.i)
    
    
    ### Now form dm/dtheta
    dg.i <- matrix(0, nrow = V + V * (V + 1)/2,
                   ncol = dim(beta.index)[1] + dim(omega.index)[1])
    
    # derivative of mean constraints
    for(k in 1:V){
      dg.i[k, which(beta.index[, 1] == k)] <- -Y[beta.index[which(beta.index[, 1] == k), 2], n]
    }
    
    for(j in 1:(V * (V+1) / 2)){
      dg.i[V + j, which(beta.index[,1] == covar.temp[j,1])] <- -Y[beta.index[which(beta.index[,1] == covar.temp[j, 1]), 2], n] * errors[covar.temp[j,2], n]
      dg.i[V + j, which(beta.index[,1] == covar.temp[j,2])] <- -Y[beta.index[which(beta.index[,1] == covar.temp[j, 2]), 2], n] * errors[covar.temp[j,1], n]
      if(Omega[covar.temp[j,1],covar.temp[j,2]] == 1){
        dg.i[V + j, dim(beta.index)[1] + which(covar.temp[j,1] == omega.index[,1] & covar.temp[j,2] == omega.index[,2])] <- -1 
      }
    }
    dg <- dg + dg.i * p[n]
  }
  
  return(solve(t(dg) %*% solve(g.var) %*% dg))
}


## Helper function for emp
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
  for(i in 1:V){
    ret <- rbind(ret, cbind(which(Omega[i:V, i]==1)+i-1, rep(i, sum(Omega[i:V, i]==1))))
  } 
  return(ret)
}



#' Maximum Likelihood Method for Linear Structural Equation Models
#' 
#' Calculate Fisher information for non-structural zeros in B and Omega matrix. The ordering of the estimates
#' in the Information matrix follows standard vec(-) convention indexing over rows first then columns. The elements of B
#' are ordered before the elements of Omega. Duplicate elements of Omega are not included.
#'
#' @param S V by V matrix corresponding to the sample covariance    
#' @param B V by V matrix with {0,1} giving structure of directed edges
#' @param Omega V by V matrix with {0,1} giving structure of bi-directed edges
#' @param B.hat V by V matrix giving estimated edges weights for directed edges
#' @param Omega.hat V by V matrix giving estimated edge weights for bi-directed edges
#' @param type string describing which Fisher Information to calculate. Options are "expected" or "observed".
#' @return The Fisher information matrix as derived by Fox and Drton 2014
#' @export
getRicfInfo <- function(S, B, Omega,  B.hat, Omega.hat, type = "expected")
{
  V <- dim(B)[1]
  ### Create P Matrix ####
  P <- matrix(0, nrow = V^2, ncol = sum(B))
  temp <- matrix(c(which(c(B)==1),1:sum(B)), nrow = sum(B), ncol = 2)
  P[temp] <- 1
  Q <- matrix(0, nrow = V^2, ncol = (sum(Omega)-V)/2 + V)
  
  ## TO-DO: Possibly use subsetting with direct insertion
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
  
  if(type == "expected"){
    I[1:sum(B), 1:sum(B)] <- t(P) %*% (S %x% omega.inv + (eye.minus.b.inv %x% t(eye.minus.b.inv)) %*% Tr) %*% P
    I[-c(1:sum(B)),-c(1:sum(B))] <- 1/2 * t(Q) %*% (omega.inv %x% omega.inv) %*% Q
    I[c(1:sum(B)),-c(1:sum(B))] <- t(P) %*% ( eye.minus.b.inv %x% omega.inv) %*% Q
    I[-c(1:sum(B)),c(1:sum(B))] <- t(Q) %*% (t(eye.minus.b.inv) %x% omega.inv) %*% P
  } else if (type == "observed"){
    back_term <- omega.inv %*% (diag(rep(1, V)) - BHat) %*% S %*% t((diag(rep(1, V)) - BHat)) %*% omega.inv
    I[1:sum(B), 1:sum(B)] <- t(P) %*% (S %x% omega.inv + (eye.minus.b.inv %x% t(eye.minus.b.inv)) %*% Tr) %*% P
    # notice difference between this and expected
    I[-c(1:sum(B)),-c(1:sum(B))] <- 1/2 * t(Q) %*% 
      ((omega.inv %x% omega.inv) -
         (omega.inv %x% back_term) - 
         (back_term %x% omega.inv)) %*% Q
    
    I[c(1:sum(B)),-c(1:sum(B))] <- t(P) %*% ((S %*% t(diag(rep(1, V)) - BHat) %*% eye.minus.b.inv) %x% omega.inv) %*% Q
    I[-c(1:sum(B)),c(1:sum(B))] <- t(I[c(1:sum(B)),-c(1:sum(B))])
  } else {
    stop("Argument type not valid")
  }
  
  return(I)
}
