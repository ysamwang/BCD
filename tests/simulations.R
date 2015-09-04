#### Simulation Studies for Paper ####
library(MASS)
set.seed(1001)
p <- 10
n <- 500
k <- 5
d <- .2
b <- d/2
B <- matrix(0, nrow = p, ncol = p)
Omega <- diag(rep(1, p))
for(i in 1:(k-1))
{
  B[i+1, i] <- 1
}
B[1, k] <- 1

for(i in 2:p)
  {
  for(j in 1:(i-1))
  {
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

reorder <- sample(p)
B <- B[reorder, reorder]
Omega <- Omega[reorder, reorder]
B.true <- matrix(rnorm(p^2), nrow = p) * B
Omega.true <- matrix(rnorm(p^2), nrow = p)
Omega.true[lower.tri(Omega.true, diag = F)] <- t(Omega.true)[lower.tri(Omega.true, diag = F)] 
Omega.true <- Omega.true * Omega
for(i in 1:p){
  Omega.true[i,i] <- sum(abs(Omega.true[i, -i])) + 1 + rchisq(1, df = 1)
}
sigma <- solve(diag(rep(1,p)) - B.true) %*% Omega.true %*% t(solve(diag(rep(1,p)) - B.true))
Y <- t(mvrnorm(n = n, mu = rep(0, p), Sigma = sigma))
Y <- Y - rowMeans(Y)

out <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
            OmegaInit = NULL, sigConv = 1, maxIter =500)
mean((out$SigmaHat - sigma)^2)
