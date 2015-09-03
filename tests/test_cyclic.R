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


out <- ricf(B = B, Omega = Omega, Y = Y, BInit = B.initial,
            OmegaInit = Omega.initial, sigConv = 0, maxIter = 1000)
sum((out$B_hat - B.true)^2)
sum((out$Omega_hat - Omega.true)^2)

sum((B.initial - B.true)^2)
sum((Omega.initial - Omega.true)^2)



sigma_hat <- solve(diag(rep(1,V)) - out$B_hat) %*% out$Omega_hat %*% t(solve(diag(rep(1,V)) - out$B_hat))
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



# P <- diag(rep(1, 5))
# P[, 1] <- P[1, ] <- c(0,0,0,1,0)
# P[, 4] <- P[4, ] <- c(1,0,0,0,0)
# O <- P %*% O %*% P
# L <- P %*% L %*% P
# X <- X %*% P
# Linit <- P %*% Linit %*% P
# Oinit <- P %*% Oinit %*% P
# testPerm <- ricfR(L = L, O = O, X  = X, Linit = LInit, Oinit = OInit, sigconv=FALSE, tol=10^(-6),
#               maxiter=1, out="none", maxkap = 10^6, B = NULL)
# 
# Omega <- P %*% Omega %*% P
# B <- P %*% B %*% P
# Y <- P %*% Y
# B.initial <- P %*% B.initial %*% P
# Omega.initial <- P %*% Omega.initial %*% P
# 
# out <- ricf(B = B, Omega = Omega, Y = Y, BInit = B.initial,
#             OmegaInit = Omega.initial, sigConv = 0, maxIter = 1)
# 
microbenchmark(ricf(B = B, Omega = Omega, Y = Y, BInit = B.initial, OmegaInit = Omega.initial))
microbenchmark(ricfR(L = L, O = O, X  = X, Linit = Linit, Oinit = Oinit, sigconv=FALSE, tol=10^(-6),
                     maxiter=500, out="none", maxkap = 10^6, B = NULL))