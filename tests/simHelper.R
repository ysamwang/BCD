#### Simulation Helper File #####
library(MASS)
library(BCD)
library(sem)
library(microbenchmark)

#### Functions for SEM ####

# returns covar argument for sem::specifyModel function
# Takes Omega Matrix of 1/0 specifying structural 0's in Bidirected Graph
getCovar <- function(Omega){
  p <- nrow(Omega)
  covar <- c()
  for(i in 1:p)
  {
    for(j in 1:i){
      if(i == j){
        covar <- c(covar, paste("x", i, sep = ""))
      } else {
      if(Omega[i,j] ==1){
        covar <- c(covar, paste(paste("x", i, sep =""), paste("x", j, sep=""), sep = ", "))
      }
    }
    
    }
  }
  return(covar)
}

# returns text argument for sem::specifyModel function
# Takes B Matrix of 1/0 specifying structural 0's in Directed Graph
# B_ij = 1 implies a directed edge from j -> i
getMod <- function(B){
  p <- nrow(B)
  mod <- ""
  first <- 1
  for(i in 1:p){
      for(j in 1:p){
        if(B[i,j] == 1){
          if(first){
            mod <- paste("x", j, "->", "x",i,", beta",i,j, sep = "")
            first <- 0 
          } else {
          mod <- paste(mod, "\n ","x", j, "->", "x",i,", beta",i,j, sep = "")
        }
      }
    }
  }
  return(mod)
}


# wrapper function for sem's specifyModel
# takes spec
specifyModelSEM <- function(B, Omega){
  return(sem::specifyModel(text = getMod(B), exog.variances = T, endog.variances = T, covs = getCovar(Omega), quiet = T))
}

#### Mixed Method ####

# returns text argument for sem::specifyModel
# takes a ricf object out, B matrix, and Omega matrix
specifyModelSEM_Mixed <- function(out, B, Omega)
{
  p <- nrow(B)
  mod <- ""
  first <- 1
  for(i in 1:p){
    for(j in 1:p){
      if(B[i,j] == 1){
        if(first){
          mod <- paste("x", j, "->", "x",i,", beta",i,j,",", out$BHat[i,j], sep = "")
          first <- 0 
        } else {
          mod <- paste(mod, "\n ","x", j, "->", "x",i,", beta",i,j,",", out$BHat[i,j], sep = "")
        }
      }
      if(Omega[i,j] == 1 & i <= j){
        if(first){
          mod <- paste("x", j, "<->", "x",i,", omega",i,j,",", out$OmegaHat[i,j], sep = "")
          first <- 0 
        } else {
          mod <- paste(mod, "\n ","x", j, "<->", "x",i,", omega",i,j,",", out$OmegaHat[i,j], sep = "")
        }
      }
    }
  }
  return(mod)
}


# returns sem object af
# Takes B, Omega and Y
mixedMethod <- function(B, Omega, Y, ricfIter = 10, ricfTol = 1e-6){
  p <- dim(B)[1]
  out <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
       OmegaInit = NULL, sigConv = 0, maxIter = ricfIter, msgs = FALSE, omegaInitScale = .9, tol = ricfTol)
  mod <- sem::specifyModel(text = specifyModelSEM_Mixed(out, B, Omega), quiet = T)
  d <- as.data.frame(t(Y))
  names(d) <- paste("x", c(1:p), sep = "")
  ret <- sem::sem(model = mod, data = d)
  return(ret)
}




##### Simulation #####
do.one <- function(p, n, k, d, b, times = 5, tol = .5){

  
  #### Generate Model ####
  
  # Setup B and Omega matrices
  B <- matrix(0, nrow = p, ncol = p)
  Omega <- diag(rep(1, p))
  
  # Generate cycle from 1 -> 2 -> ... K -> 1
  if(k > 0){
    for(i in 1:(k-1))
    {
      B[i+1, i] <- 1
    }
    B[1, k] <- 1
  }

  # fill in remaining edges  
  for(i in 2:p)
  {
    for(j in 1:(i-1))
    {
      if(!B[j,i]){
      U <- runif(1)
      if(U < d){
        B[j, i] <-1
      } else {
        if(U < b + d)
          { 
          Omega[i,j] <- Omega[j,i] <- 1  
          }
        }
      }
    }
  }
  
  # reorder the vertices
  reorder <- sample(p)
  B <- B[reorder, reorder]
  Omega <- Omega[reorder, reorder]

  # Sample edge weights as given in the paper
  B.true <- matrix(rnorm(p^2), nrow = p) * B
  Omega.true <- matrix(rnorm(p^2), nrow = p)
  Omega.true[lower.tri(Omega.true, diag = F)] <- t(Omega.true)[lower.tri(Omega.true, diag = F)] 
  Omega.true <- Omega.true * Omega
  
  # ensure Omega.true is PD by making it diagonally dominant
  for(i in 1:p){
    Omega.true[i,i] <- sum(abs(Omega.true[i, -i])) + 1 + rchisq(1, df = 1)
  }
  sigma <- solve(diag(rep(1,p)) - B.true) %*% Omega.true %*% t(solve(diag(rep(1,p)) - B.true))
  
  # Sample data from multivariate normal and make mean 0
  Y <- t(MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = sigma))
  Y <- Y - rowMeans(Y)

  


   #### Run sem ####
  dat <- as.data.frame(t(Y))
  names(dat) <- paste("x", c(1:p), sep = "")
  mod <- specifyModelSEM(B, Omega)
  time.sem <- try(microbenchmark::microbenchmark(out.sem <- sem(model = mod, data = dat), times = times),
                  silent = T)
  
  #### Run RICF ####
  time.ricf <- try(microbenchmark::microbenchmark(out.ricf <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                                                                   OmegaInit = NULL, sigConv = 0, maxIter = 5000,
                                                                   msgs = FALSE, omegaInitScale = .9)
                                                  , times = times), silent = T)
  if(class(time.ricf)[1] == "try-error")
  {
    time.ricf <- list(time = 0)
    out.ricf <- list(Converged = 0)
  }
  
  if(class(time.sem)[1] == "try-error")
  {
    time.sem <- list(time = 0)
    out.sem <- list(convergence = 0)
  }  else if(any(eigen(out.sem$C)$values < 0 )){
    time.sem <- list(time = 0)
    out.sem <- list(convergence = 0)
  }
  
  
   #### Run Mixed ####
  time.mixed <- try(microbenchmark::microbenchmark(out.mixed <- mixedMethod(B, Omega, Y), times = times), silent = T)
  
  if(class(time.mixed)[1] == "try-error")
  {
    time.mixed <- list(time = 0)
    out.mixed <- list(convergence = 0)
  } else if(any(eigen(out.mixed$C)$values < 0 )){
    time.mixed <- list(time = 0)
    out.mixed <- list(convergence = 0)
  }
  
  #### Post Processing
  
  if(out.ricf$Converged & out.sem$convergence){
    err <- (sum(abs(out.ricf$BHat - out.sem$A)) + sum(abs(out.ricf$OmegaHat - out.sem$P ))) / sum(B + Omega)    
    agree <- err < tol
  } else {
    err <- -1
    agree <- 0
  }
  if(is.na(agree))
  {}

  ret <- list(ricfTime = mean(time.ricf$time), ricfConv = out.ricf$Converged,
                semTime = mean(time.sem$time), semConv = out.sem$convergence,
                mixedTime = mean(time.mixed$time), mixedConv = out.mixed$convergence,
                agree = agree, B = B, Omega = Omega, B.true = B.true, Omega.true = Omega.true, Y = Y)
  
  return(ret)
}


#### Generate Model ####
generateModel <-function(v, n, k, d, b){
      # Setup B and Omega matrices
      B <- matrix(0, nrow = v, ncol = v)
      Omega <- diag(rep(1, v))
      
      # Generate cycle from 1 -> 2 -> ... K -> 1
      if(k > 0){
        for(i in 1:(k-1))
        {
          B[i+1, i] <- 1
        }
        B[1, k] <- 1
      }
      
      # fill in remaining edges  
      for(i in 2:v)
      {
        for(j in 1:(i-1))
        {
          if(B[j,i]!= 1){
            U <- runif(1)
            if(U < d){
              B[j, i] <-1
            } else {
              if((B[i,j]!= 1) & (U < b + d))
              { 
                Omega[i,j] <- Omega[j,i] <- 1  
              }
            }
          }
        }
      }
      
      # reorder the vertices
      reorder <- sample(v)
      B <- B[reorder, reorder]
      Omega <- Omega[reorder, reorder]
      
      # Sample edge weights as given in the paper
      B.true <- matrix(rnorm(v^2), nrow = v) * B
      Omega.true <- matrix(rnorm(v^2), nrow = v)
      Omega.true[lower.tri(Omega.true, diag = F)] <- t(Omega.true)[lower.tri(Omega.true, diag = F)] 
      Omega.true <- Omega.true * Omega
      
      # ensure Omega.true is PD by making it diagonally dominant
      for(i in 1:v){
        Omega.true[i,i] <- sum(abs(Omega.true[i, -i])) + 1 + rchisq(1, df = 1)
      }
      sigma <- solve(diag(rep(1,v)) - B.true) %*% Omega.true %*% t(solve(diag(rep(1,v)) - B.true))
      
      # Sample data from multivariate normal and make mean 0
      Y <- t(MASS::mvrnorm(n = n, mu = rep(0, v), Sigma = sigma))
      Y <- Y - rowMeans(Y)
      return(list(Y= Y, B = B, Omega = Omega, B.true = B.true, Omega.true = Omega.true))
}

