#### BCD Test File ####

# cyclic B
library(MASS)
library(microbenchmark)
set.seed(1001)
n <- 20
V <- 5
B <- matrix(c(0,1,1,0,0,
              0,0,1,0,0,
              0,0,0,1,0,
              0,0,0,0,1,
              1,0,1,0,0), nrow = V, byrow = T)

Omega <- diag(rep(1,V))
Omega[1,3] <- Omega[3,1] <- 1

B.true <- matrix(rnorm(V^2), nrow = V) * B
Omega.true <- ifelse(Omega,runif(V^2, -1, 1),0) + diag(rep(4, V)) 
sigma = solve(diag(rep(1, V)) - B.true) %*% Omega.true %*% t(solve(diag(rep(1, V)) - B.true))

Y <- t(mvrnorm(n = n, mu = rep(0, V), Sigma = sigma))
Y <- Y - rowMeans(Y)
B.initial <- matrix(rnorm(V^2), nrow = V) * B

Omega.initial <- matrix(rnorm(V^2), nrow = V)
Omega.initial <- Omega.initial %*% t(Omega.initial)
Omega.initial <- Omega.initial * Omega
all(eigen(Omega.initial)$values > 0)
sigma.initial = solve(diag(rep(1, V)) - B.initial) %*% Omega.initial %*% t(solve(diag(rep(1, V)) - B.initial))
all(eigen(sigma.initial)$values > 0)


sum((B.initial - B.true)^2)
sum((Omega.initial - Omega.true)^2)


out <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
            OmegaInit = NULL, sigConv = 0, maxIter = 5000)
sum((out$BHat - B.true)^2)
sum((out$OmegaHat - Omega.true)^2)
out$Iter

out <- ricf(B = B, Omega = Omega, Y = Y, BInit = B.initial,
            OmegaInit = Omega.initial, sigConv = 0, maxIter = 5000)
sum((out$BHat - B.true)^2)
sum((out$OmegaHat - Omega.true)^2)
out$Iter


sum((B.initial - B.true)^2)
sum((Omega.initial - Omega.true)^2)

microbenchmark(ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                                   OmegaInit = NULL, sigConv = 0, maxIter = 5000))
microbenchmark(ricf(B = B, Omega = Omega, Y = Y, BInit = B.initial,
                    OmegaInit = Omega.initial, sigConv = 0, maxIter = 5000))


sigma_hat <- solve(diag(rep(1,V)) - out$BHat) %*% out$OmegaHat %*% t(solve(diag(rep(1,V)) - out$BHat))
sigma.true <- solve(diag(rep(1,V)) - B.true) %*% Omega.true %*% t(solve(diag(rep(1,V)) - B.true))



L <- t(B)
O <- Omega
X <- t(Y)
Linit <- t(B.initial)
Oinit <- Omega.initial

test <- ricfR(L = L, O = O, X  = X, Linit = Linit, Oinit = Oinit, sigconv=FALSE, tol=10^(-6),
                  maxiter=500, out="none", maxkap = 10^6, B = NULL)

sum((out$B_hat + test$Bhat)^2)
sum((out$Omega_hat - test$Omegahat)^2)


microbenchmark(ricf(B = B, Omega = Omega, Y = Y, BInit = B.initial, OmegaInit = Omega.initial))
microbenchmark(ricfR(L = L, O = O, X  = X, Linit = Linit, Oinit = Oinit, sigconv=FALSE, tol=10^(-6),
                     maxiter=500, out="none", maxkap = 10^6, B = NULL))