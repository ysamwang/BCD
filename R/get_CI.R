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
#' @param type the type of variances to return. options are "asymptotic" for the asymptotic
#' variance as defined in Qin and Lawless 1994 or "profile" for profile confiendence intervals
#' @return The (scaled by n) estimated covariance matrix or profile confidence intervals
#' @export
var.el <- function(Y, B, Omega, B.hat, Omega.hat = NULL, p = NULL, cutoff = NULL, grid = 5e-3, type = "asymptotic")
{
  if(type == "asymptotic"){
    # Index of non-zero B's in proper order (column by column)
    beta.index <- .getBetaIndex(B)
    # Index of non-zero Omega's
    omega.index <- .getOmegaIndex(Omega)
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
    
    return(solve(t(dg) %*% solve(g.var) %*% dg) /n)
    
  } else if(type == "profile-fixed"){
    
    columns <- matrix(c(1:V), nrow = V, ncol = V, byrow = T)
    rows <- matrix(c(1:V), nrow = V, ncol = V)
    beta.indices <- cbind(rows[which(B==1)], columns[which(B==1)])
    
    ret <- matrix(0, nrow = dim(beta.indices)[1], ncol = 5)
    ret[, 2] <- B.hat[which(B == 1)]
    
    for(i in 1:dim(beta.indices)[1]){
      ret[i, -2] <- .profile.like.fixed(Y, B, Omega, B.hat, cutoff, index = i, grid = grid)
    }
    
    colnames(ret) <- c("Lower", "Estimate", "Upper", "Lower-LRT", "Upper-LRT")
    return(ret)
  }  else if(type == "profile-full"){
    columns <- matrix(c(1:V), nrow = V, ncol = V, byrow = T)
    rows <- matrix(c(1:V), nrow = V, ncol = V)
    beta.indices <- cbind(rows[which(B==1)], columns[which(B==1)])
    
    ret <- matrix(0, nrow = dim(beta.indices)[1], ncol = 5)
    ret[, 2] <- B.hat[which(B == 1)]
    method = "Nelder"
    
    for(i in 1:dim(beta.indices)[1]){
      ret[i, -2] <- .profile.like.full(Y, B, Omega, B.hat, cutoff, grid, beta.indices[i,1], beta.indices[i,2], method = method)
    }
    
    colnames(ret) <- c("Lower", "Estimate", "Upper", "Lower-LRT", "Upper-LRT")
    return(ret)
  } else {
    stop("Invalid type given. Options are: asymptotic, profile-fixed, profile-full")
  }
  
}


## Helper function for empirical likelihood variance
.getBetaIndex <- function(B)
{
  ret <- NULL
  for(i in 1:dim(B)[1]){
      ret <- rbind(ret, cbind(which(B[,i]==1) , rep(i, sum(B[,i]))))
  }
  return(ret)
}

## Helper function for empirical likelihood variance
.getOmegaIndex <- function(Omega){ 
  ret <- NULL
  V <- dim(Omega)[1]
  for(i in 1:V){
    ret <- rbind(ret, cbind(which(Omega[i:V, i]==1)+i-1, rep(i, sum(Omega[i:V, i]==1))))
  } 
  return(ret)
}


.profile.like.fixed <- function(Y, B, Omega, B.hat, cutoff, index, grid = 1e-3) {
  tol = 1e-6
  maxInnerIter = 100
  mu.hat = NULL
  num_dual_vars = V + sum(Omega == 0)/2
  val <- B.hat[which(B==1)]
  point.est <- val[index]
  
  fitted_mod <- sem_el_fit_weights(val, y_r = Y, omega_r = Omega, b_r = B,
                                   dual_r = rep(0, num_dual_vars), tol = tol, max_iter = maxInnerIter, meanEst = !is.null(mu.hat))
  
  
  # get upper bound
  while(-2 * sum(log(1/fitted_mod$d * n) ) < cutoff){
    val[index] <- val[index] + grid
    fitted_mod <- sem_el_fit_weights(val, y_r = Y, omega_r = Omega, b_r = B,
                                     dual_r = rep(0, num_dual_vars), tol = tol, max_iter = maxInnerIter, meanEst = !is.null(mu.hat))  
  }
  upper <- val[index]- grid
  upper.chisq <- -2 * sum(log(1/fitted_mod$d * n) )
  
  # get lower bound
  val[index] <- point.est
  fitted_mod <- sem_el_fit_weights(val, y_r = Y, omega_r = Omega, b_r = B,
                                   dual_r = rep(0, num_dual_vars), tol = tol, max_iter = maxInnerIter, meanEst = !is.null(mu.hat))
  
  while(-2 * sum(log(1/fitted_mod$d * n))  < cutoff){
    val[index] <- val[index] - grid
    fitted_mod <- sem_el_fit_weights(val, y_r = Y, omega_r = Omega, b_r = B,
                                     dual_r = rep(0, num_dual_vars), tol = tol, max_iter = maxInnerIter, meanEst = !is.null(mu.hat))  
  }
  lower <- val[index] + grid
  lower.chisq <- -2 * sum(log(1/fitted_mod$d * n) )
  
  return(c(lower, upper, upper.chisq, lower.chisq))
}




.profile.like.full <- function(Y, B, Omega, B.hat, cutoff, grid, row.ind, col.ind, method = method) {
  tol = 1e-6
  maxInnerIter = 100
  meanEst = 0
  
  
  B.mod <- B 
  B.mod[row.ind, col.ind] <- 0
  
  init_val <- B.hat[B.mod == 1]
  
  # record minimizer point
  point.est <- B.hat[row.ind, col.ind]
  
  
  # Lower bound
  LRT <- 0
  while(LRT < cutoff){
    
    B.hat[row.ind, col.ind] <- B.hat[row.ind, col.ind] - grid  
    optim_out <- optim(par = init_val, fn = sem_el_fit_obj_one_fixed, gr = NULL,
                       y_r = Y, omega_r = Omega, b_r = B.mod,
                       dual_r = rep(0, V + sum(Omega == 0)/2),
                       tol = tol, max_iter = maxInnerIter,
                       meanEst = meanEst, b_fixed = B.hat[row.ind, col.ind], row_ind = row.ind, col_ind = col.ind,
                       method = method, control = list(fnscale = -1))
    
    B.hat[B.mod == 1] <- optim_out$par
    fitted_mod <- sem_el_fit_weights(c(B.hat[B == 1]), y_r = Y, omega_r = Omega, b_r = B,
                                     dual_r = rep(0, V + sum(Omega == 0)/2),
                                     tol = tol, max_iter = maxInnerIter, meanEst = meanEst)
    LRT <- -2 * sum(log(n / c(fitted_mod$d)))
  }
  lower <- B.hat[row.ind, col.ind] + grid
  lower.lrt <- LRT
  
  
  
  # Upper bound
  LRT <- 0
  B.hat[row.ind, col.ind] <- point.est
  
  while(LRT < cutoff){
    
    B.hat[row.ind, col.ind] <- B.hat[row.ind, col.ind] + grid  
    optim_out <- optim(par = init_val, fn = sem_el_fit_obj_one_fixed, gr = NULL,
                       y_r = Y, omega_r = Omega, b_r = B.mod,
                       dual_r = rep(0, V + sum(Omega == 0)/2),
                       tol = tol, max_iter = maxInnerIter,
                       meanEst = meanEst, b_fixed = B.hat[row.ind, col.ind], row_ind = row.ind, col_ind = col.ind,
                       method = method, control = list(fnscale = -1))
    
    B.hat[B.mod == 1] <- optim_out$par
    fitted_mod <- sem_el_fit_weights(c(B.hat[B == 1]), y_r = Y, omega_r = Omega, b_r = B,
                                     dual_r = rep(0, V + sum(Omega == 0)/2),
                                     tol = tol, max_iter = maxInnerIter, meanEst = meanEst)
    LRT <- -2 * sum(log(n / c(fitted_mod$d)))
  }
  upper <- B.hat[row.ind, col.ind] - grid
  upper.lrt <- LRT
  return(c(lower, upper, lower.lrt, upper.lrt))
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
#' @return The inverse (scaled by n) Fisher information matrix as derived by Fox and Drton 2014 or the Huber-White misspecified model covariance estimate
#' @export
var.ricf <- function(Y, B, Omega,  B.hat, Omega.hat, type = "expected")
{
  
  # Model based SE's with either expected or observed information
  if(type == "expected" | type == "observed"){
    V <- dim(B)[1]
    I <- matrix(0, nrow = sum(B) + (sum(Omega)-V)/2 + V,
                ncol = sum(B) + (sum(Omega)-V)/2 + V)
    n <- dim(Y)[2]
    S <- Y %*% t(Y) / n
    
  
    #### Setup permutation matrix P for B
    P <- matrix(0, nrow = V^2, ncol = sum(B))
    temp <- matrix(c(which(c(B)==1),1:sum(B)), nrow = sum(B), ncol = 2)
    P[temp] <- 1
    
    
    ### Setup permutation matrix Q for Omega
    Omega.temp <- matrix(0, nrow = V, ncol = V)
    Omega.temp[lower.tri(Omega.temp, diag = T) & Omega == 1] <- c(1:((sum(Omega) - V)/2 + V))
    Omega.temp <- Omega.temp + t(Omega.temp) - diag(diag(Omega.temp))
    Q <- matrix(0, nrow = V^2, ncol = (sum(Omega)-V)/2 + V)
    Q[cbind(c(1:V^2), c(Omega.temp))] <- 1
    
    
    ### Setup Transpositon matrix such that vec(t(I-B)) = Tr %*% vec(I-B)
    Tr <- matrix(0, nrow = V^2, ncol = V ^2)
    Tr[matrix(c(1:V^2, (rep(c(1:V), V)-1) *V + rep(c(1:V), each = V)), ncol = 2)] <- 1
    
    ### precompute quantities of interest
    omega.inv <- solve(Omega.hat)
    eye.minus.b.inv <- solve(diag(rep(1, V)) - B.hat)
    
    
    if(type == "expected"){
      ## Expected Information where we assume E(S) = Sigma
      
      # information for beta elements
      I[1:sum(B), 1:sum(B)] <- t(P) %*% (S %x% omega.inv + (eye.minus.b.inv %x% t(eye.minus.b.inv)) %*% Tr) %*% P
      
      # information for omega elements
      I[-c(1:sum(B)),-c(1:sum(B))] <- 1/2 * t(Q) %*% (omega.inv %x% omega.inv) %*% Q
      
      # co-information
      I[c(1:sum(B)),-c(1:sum(B))] <- t(P) %*% ( eye.minus.b.inv %x% omega.inv) %*% Q
      I[-c(1:sum(B)),c(1:sum(B))] <- t(I[c(1:sum(B)),-c(1:sum(B))])
      
    } else if (type == "observed"){
      ## Observed Information where we do not assume E(S) = Sigma
      
      # precompute term
      back_term <- omega.inv %*% (diag(rep(1, V)) - B.hat) %*% S %*% t((diag(rep(1, V)) - B.hat)) %*% omega.inv
      
      # information for beta elements
      I[1:sum(B), 1:sum(B)] <- t(P) %*% (S %x% omega.inv + (eye.minus.b.inv %x% t(eye.minus.b.inv)) %*% Tr) %*% P
      
      # information for omega elements (notice difference between this and expected info above)
      I[-c(1:sum(B)),-c(1:sum(B))] <- -1/2 * t(Q) %*%((omega.inv %x% omega.inv) - (omega.inv %x% back_term) - 
           (back_term %x% omega.inv)) %*% Q
      
      # co-information for omega elements
      I[c(1:sum(B)),-c(1:sum(B))] <- t(P) %*% ((S %*% t(diag(rep(1, V)) - B.hat) %*% omega.inv) %x% omega.inv) %*% Q
      I[-c(1:sum(B)),c(1:sum(B))] <- t(I[c(1:sum(B)),-c(1:sum(B))])
    }
    
    # return inverse of FI
    return(solve(I)/n)
    
  } else if (type == "sandwich") {
    #### Sandwich Variance based on misspecified model ####
    
    V <- dim(B)[1]
    
    S <- Y %*% t(Y) / dim(Y)[2]
    
    ## Setup permutation matrix P for B
    P <- matrix(0, nrow = V^2, ncol = sum(B))
    temp <- matrix(c(which(c(B)==1),1:sum(B)), nrow = sum(B), ncol = 2)
    P[temp] <- 1
    
    
    ## Setup permutation matrix Q for Omega
    Omega.temp <- matrix(0, nrow = V, ncol = V)
    Omega.temp[lower.tri(Omega.temp, diag = T) & Omega == 1] <- c(1:((sum(Omega) - V)/2 + V))
    Omega.temp <- Omega.temp + t(Omega.temp) - diag(diag(Omega.temp))
    Q <- matrix(0, nrow = V^2, ncol = (sum(Omega)-V)/2 + V)
    Q[cbind(c(1:V^2), c(Omega.temp))] <- 1
    
    
    ## Setup Transpositon matrix such that vec(t(I-B)) = Tr %*% vec(I-B)
    Tr <- matrix(0, nrow = V^2, ncol = V ^2)
    Tr[matrix(c(1:V^2, (rep(c(1:V), V)-1) *V + rep(c(1:V), each = V)), ncol = 2)] <- 1
    
    ## precompute quantities of interest
    omega.inv <- solve(Omega.hat)
    eye.minus.b.inv <- solve(diag(rep(1, V)) - B.hat)
    omega.inv.eye.minus.b <- omega.inv %*% (diag(rep(1, V)) - B.hat)
    
    K <- J <-  matrix(0, nrow = sum(B) + (sum(Omega)-V)/2 + V,
                      ncol = sum(B) + (sum(Omega)-V)/2 + V)
    
    for(i in 1:dim(Y)[2]){
      l <- c(t(P) %*% c(omega.inv.eye.minus.b %*% (Y[,i] %*% t(Y[,i])) - t(eye.minus.b.inv)),
             t(Q) %*% c(omega.inv - omega.inv.eye.minus.b %*% Y[,i] %*% t(Y[,i]) %*% t(omega.inv.eye.minus.b)) / (-2) ) 
      K <- K + l %*% t(l)
    }
    K <- K / dim(Y)[2]
    
    # precompute term
    back_term <- omega.inv %*% (diag(rep(1, V)) - B.hat) %*% S %*% t((diag(rep(1, V)) - B.hat)) %*% omega.inv
    
    # information for beta elements
    J[1:sum(B), 1:sum(B)] <- t(P) %*% (S %x% omega.inv + (eye.minus.b.inv %x% t(eye.minus.b.inv)) %*% Tr) %*% P
    
    # information for omega elements (notice difference between this and expected info above)
    J[-c(1:sum(B)),-c(1:sum(B))] <- -1/2 * t(Q) %*%((omega.inv %x% omega.inv) - (omega.inv %x% back_term) - 
                                                      (back_term %x% omega.inv)) %*% Q
    
    # co-information for omega elements
    J[c(1:sum(B)),-c(1:sum(B))] <- t(P) %*% ((S %*% t(diag(rep(1, V)) - B.hat) %*% omega.inv) %x% omega.inv) %*% Q
    J[-c(1:sum(B)),c(1:sum(B))] <- t(J[c(1:sum(B)),-c(1:sum(B))])
    
    # return estimate of variance
    return(solve(J) %*%K %*% t(solve(J)) / n)
  } else {
    stop("Invalid type given. Options are: observed, expected, sandwich")
  }

}


