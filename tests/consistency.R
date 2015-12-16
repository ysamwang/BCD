### Simulation Goals 
#
# Question: 
# Under error misspecification, are the estimates still consistent and what does
# finite sample bias look like
#
# Approach:
# Simulate errors from lognormal distribution and then examine bias
# and MSE as sample size grows
#
#

library('BCD')
library('MASS')

set.seed(123)
V <- 8
n <- 1000
k <- 3
b <- .3
d <- .3
dist <- "lognormal"
mod <- generateModel(v = V, n = n, k = k, b = b, d = d, errorDist = dist)

Sigma <- solve(diag(rep(1, V)) - mod$B.true) %*% mod$Omega.true %*% t(solve(diag(rep(1, V)) - mod$B.true))

mult.factor <- solve(diag(rep(1,V)) - mod$B.true) %*% t(chol(mod$Omega.true))

omega.true <- mod$Omega.true[lower.tri(mod$Omega.true, diag = T)][which(mod$Omega[lower.tri(mod$Omega, diag = T)]==1)]
b.true <- mod$B.true[which(mod$B==1)]
sim.size <- 500
n.list <- c(100, 200,500, 1000, 5000, 10000)

for(type in c("lognormal_bias")){ 
  for(n.index in 1:length(n.list)){
    b.mse.like <-  o.mse.like <- b.bias.like <- o.bias.like <- rep(0, sim.size)
    bias.record <- mse.record <- matrix(0, nrow = length(n.list), ncol = sum(mod$B) + (sum(mod$Omega) - V)/2 + V)
    n <- n.list[n.index]
    cat(paste("\n===",n, type,"\n"))
    for(i in 1:sim.size){
      cat(paste(i,": "))
      
      
      if(type == "gauss"){
        Y <- t(MASS::mvrnorm(n = n, Sigma = Sigma, mu = rep(0, V)))
        Y <- Y - rowMeans(Y)
      } else {
        mu <- 0
        sigma <- 1 
        eps <- matrix(exp(rnorm(n * V, mean = mu, sd = sigma)),
                      nrow = V, ncol = n) - exp(mu + sigma^2/2)
        eps <- eps / sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2))
        Y <- mult.factor %*% eps
      }
      
      
      
      
      #### RICF Procedure ####
      out_ricf <- ricf(B = mod$B, Omega = mod$Omega, Y = Y, BInit = NULL,
                       OmegaInit = NULL, sigConv = 0, maxIter = 10000,
                       msgs = FALSE, omegaInitScale = .9, tol = 1e-6)
      
      
      beta.like <- out_ricf$BHat[which(mod$B==1)]
      omega.like <- out_ricf$OmegaHat[lower.tri(mod$Omega.true, diag = T) & mod$Omega ==1]
      
      b.mse.like[i,] <-  (beta.like - b.true)^2
      b.bias.like[i,] <- (beta.like - b.true)
      o.mse.like[i,] <-  (omega.like - omega.true)^2
      o.bias.like[i,] <- (omega.like - omega.true)
      
    }
    bias.record[n.index, ] <- c(b.bias.like, o.bias.like)/ sim.size
    mse.record[n.index, ] <- c(b.mse.like, o.mse.like)/ sim.size
  }
}
