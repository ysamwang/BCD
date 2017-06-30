#' Generate Data from Structural Equation Model
#' 
#' 
#' Samples a random Graph and Data
#' 
#'  
#' @param v number of nodes
#' @param n number of samples
#' @param k the size of the largest cycle
#' @param d the   
#' @return \item{sigmaHat}{estimated covariance matrix at convergence}
#'    \item{bHat}{estimated B matrix (edge weights for directed edges) at converegence}
#'    \item{omegaHat}{estimated Omega (edge weights for bi-directed edges) at convergence}
#'    \item{iterations}{number of iterations until convergence. a single iteration is considered
#'    a pass through all nodes}
#'    \item{converged}{boolean on whether or not the algorithm converged before the max iterations}
generateModel <-function(v, n, k, d, b, errorDist = "gauss"){
  # Setup B and Omega matrices
  
  
  B <- matrix(0, nrow = v, ncol = v)
  Omega <- diag(rep(1, v))
  
  # Generate cycle from 1 -> 2 -> ... K -> 1
  if(k > 1){
    for(i in 1:(k-1))
    {
      B[i+1, i] <- 1
    }
    B[1, k] <- 1
  }
  
  
  # fill in remaining edges  
  for(j in 1:(v-1))
  {
    for(i in (j+1):v)
    {
      if((B[i,j] != 1) & (B[j,i]!= 1)){
        U <- runif(1)
        if(U < d){
          B[i, j] <- 1
        } else if(U < (b + d))
        { 
          Omega[i,j] <- Omega[j,i] <- 1  
        }
      }
    }
  }
  
  # reorder the vertices
  reorder <- sample(V)
  B <- B[reorder, reorder]
  Omega <- Omega[reorder, reorder]
  
  # Sample edge weights as given in the paper
  B.true <- matrix(runif(v^2,-1,1), nrow = v) * B
  Omega.true <- matrix(runif(v^2,-1,1), nrow = v)
  Omega.true[lower.tri(Omega.true, diag = F)] <- t(Omega.true)[lower.tri(Omega.true, diag = F)] 
  Omega.true <- Omega.true * Omega
  
  # ensure Omega.true is PD by making it diagonally dominant
  for(i in 1:v){
    Omega.true[i,i] <- sum(abs(Omega.true[i, -i])) + 1 + rchisq(1, df = 1)
  }
  
  temp <- solve(diag(rep(1,V)) - B.true)
  sigma <- temp %*% Omega.true %*% t(temp)
  
  # Sample data from multivariate normal and make mean 0
  if(errorDist == "gamma")
  {
    errors <- matrix(rgamma(v * n, shape = 1, scale = 1), nrow = v, ncol = n)
    Y <- solve(diag(rep(1, v)) - B.true, t(chol(Omega.true)) %*% errors)
  }else if(errorDist == "poisson")
  {
    pois.mu <- 3
    errors <- matrix(rpois(v * n, lambda = pois.mu), nrow = v, ncol = n)
    errors <- (errors - pois.mu)/sqrt(pois.mu)
    Y <- solve(diag(rep(1, v)) - B.true, (chol(Omega.true)) %*% errors)
  } else if(errorDist == "t"){
    t.df <- 5
    errors <- matrix(rt(v * n, df = t.df) * (t.df - 2)/t.df, nrow = v, ncol = n)
    Y <- solve(diag(rep(1, v)) - B.true, t(chol(Omega.true)) %*% errors)
  } else if(errorDist == "lognormal"){
    ln.var <- 1
    ln.mu <- 0
    errors <- matrix(rlnorm(v * n, meanlog = ln.mu, sdlog = sqrt(ln.var)), nrow = v, ncol = n)
    errors <- (errors - exp(ln.mu + ln.var/2)) / sqrt((exp(ln.var) - 1) * exp(2*ln.mu + ln.var))
    Y <- solve(diag(rep(1, v)) - B.true, t(chol(Omega.true)) %*% errors)
  } else if(errorDist == "gauss") {
    Y <- t(MASS::mvrnorm(n = n, mu = rep(0, v), Sigma = sigma))
  }
  
  return(list(Y= Y, B = B, Omega = Omega, B.true = B.true,
              Omega.true = Omega.true, Sigma = sigma))
}