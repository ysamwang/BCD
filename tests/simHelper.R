#### Simulation Helper File #####

#### Functions for SEM ####

# specify covar argument given an Omega matrix
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

# specify a model argument given a B matrix
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
specifyModelSEM <- function(B, Omega){
  return(sem::specifyModel(text = getMod(B), exog.var = T, endog.var = T, covs = getCovar(Omega), quiet = T))
}


##### Simulation #####
p <- 10
n <- 500
k <- 5
d <- .2
b <- d/2

do.one <- function(p, n, k, d, b, times = 5, tol = 1e-2){

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
  B <- B[reorder, ]
  B <- B[, reorder]
  Omega <- Omega[reorder, ]
  Omega <- Omega[, reorder]
  
  # Sample edge weights as given in the paper
  B.true <- matrix(rnorm(p^2), nrow = p) * B
  Omega.true <- matrix(rnorm(p^2), nrow = p)
  Omega.true[lower.tri(Omega.true, diag = F)] <- t(Omega.true)[lower.tri(Omega.true, diag = F)] 
  Omega.true <- Omega.true * Omega
  
  for(i in 1:p){
    Omega.true[i,i] <- sum(abs(Omega.true[i, -i])) + 1 + rchisq(1, df = 1)
  }
  sigma <- solve(diag(rep(1,p)) - B.true) %*% Omega.true %*% t(solve(diag(rep(1,p)) - B.true))
  
  # Sample data from multivariate normal and make mean 0
  Y <- t(MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = sigma))
  Y <- Y - rowMeans(Y)
  
  
  ricfRres <- ricfR(O = Omega, X = t(Y), Linit = NULL, Oinit = NULL, sigconv=FALSE, B = B)
  
  time.ricf <- microbenchmark::microbenchmark(out.ricf <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                                               OmegaInit = NULL, sigConv = 0, maxIter =5000, msgs = FALSE, omegaInitScale = .9), times = times)
  if(!out.ricf$Converged){
    browser()
  }
  d <- as.data.frame(t(Y))
  names(d) <- paste("x", c(1:p), sep = "")
  
  # In case semResult throws error
  semResult <- tryCatch({time.sem <- microbenchmark::microbenchmark(out.sem <- sem(model = specifyModelSEM(B, Omega), data = d), times = times)
  1}
           , error = function(err)
             {print("sem ERROR")
             return(0)
           }, finally ={
               1
             }
           )
  
  
  if(semResult){
    err <- norm(out.ricf$BHat - out.sem$A, type = "F") + norm(out.ricf$OmegaHat - out.sem$P, type = "F")
    agree <- err < tol
    
    ret <- list(ricfTime = mean(time.ricf$time), ricfConv = out.ricf$Converged,
                semTime = mean(time.sem$time), semConv = out.sem$convergence,
                agree = err, B = B, Omega = Omega, B.true = B.true, Omega.true = Omega.true, Y = Y)
    
  } else {
    ret <- list(ricfTime = mean(time.ricf$time), ricfConv = out.ricf$Converged,
                semTime = -1, semConv = 0,
                agree = 0, B = B, Omega = Omega, B.true = B.true, Omega.true = Omega.true, Y = Y)
  }
  
  
  return(ret)
}



