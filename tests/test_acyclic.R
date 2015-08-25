#### BCD Test File ####

# True Values acyclic B
library(MASS)
library(microbenchmark)
set.seed(888)
V <- 5
B <- matrix(c(0,1,0,0,0,
              0,0,1,0,0,
              0,0,0,1,0,
              0,0,0,0,1,
              0,0,0,0,0), nrow = V, byrow = T)

Omega <- diag(rep(1,V))
Omega[1,3] <- Omega[3,1] <- 1

B.true <- matrix(rnorm(V^2), nrow = V) * B
Omega.true <- ifelse(Omega,runif(V^2, -1, 1),0) + diag(rep(4, V)) 
sigma = solve(diag(rep(1, V)) - B.true) %*% Omega.true %*% t(solve(diag(rep(1, V)) - B.true))
prod(eigen(sigma)$values > 0)

Y <- t(mvrnorm(n = 10, mu = rep(0, V), Sigma = sigma))
B.initial <- matrix(rnorm(V^2), nrow = V) * B
Omega.initial <- diag(rnorm(V))
Omega.initial[1,3] <- Omega.initial[3,1] <- min(abs(Omega.initial[1,1]), abs(Omega.initial[3,3]))/2

sum((B.initial - B.true)^2)
sum((Omega.initial - Omega.true)^2)


out <- ricf(B = B, Omega = Omega, Y = Y, BInit = B.initial,
     OmegaInit = Omega.initial, sigConv = 0, maxIter = 500)
sum((out$B_hat - B.true)^2)
sum((out$Omega_hat - Omega.true)^2)

sum((B.initial - B.true)^2)
sum((Omega.initial - Omega.true)^2)

microbenchmark(ricf(B = B, Omega = Omega, Y = Y, BInit = B.initial, OmegaInit = Omega.initial))
