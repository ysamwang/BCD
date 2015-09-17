###### Try to find something outside the model #####
library(MASS)
library(BCD)
setwd("C:/Users/ysamwang/Dropbox/Cyclic Graphs/testCases/")
V <- 5
n <- 10
counter <- 0
tries <- 0
d <- .3
b <- d/2
set.seed(1001)
while(0 < 1){
  m <- matrix(rnorm(V^2), nrow = V)
  sigma <- m %*% t(m)
  Y <- t(MASS::mvrnorm(n = n, mu = rep(0, V), Sigma = sigma))
  B <- matrix(0, V, V)
  Omega <- diag(rep(1,V))
  for(i in 2:V){
    for(j in 1:(i-1)){
      u <- runif(1)
        if(u < d)
        {
          B[i, j] <- 1
        } else if ( u < (b + d))
        {
          Omega[i,j] <- Omega[j, i] <- 1
        }
      }
    }
  out <- ricf(B = B, Omega = Omega, maxIter = 200000, Y = Y)
  tries <- tries + 1
  if(tries %% 100 == 0){
    cat(".")
  }
  if(!out$Converged){
    if(max(abs(out$BHat)) > 50000 | max(abs(out$OmegaHat))  > 50000){
      counter <- counter + 1
      cat(paste("\n Non-Convergence:", counter, "\n"))
      saveRDS(list(out2, sigma, Y, B, Omega), paste("nonConverged", counter,".RDS", sep = ""))
    }
  }
}

failure <- readRDS("nonConverged1.RDS")

out <- ricf(B = failure[[4]], Omega = failure[[5]],
            Y = failure[[3]], maxIter = 1000000)
out$Converged
out$OmegaHat

