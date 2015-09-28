#### BCD Test File ####

# True Values acyclic B
library(MASS)
library(microbenchmark)
set.seed(100)
V <- 5
n <- 10
B <- matrix(c(0,1,0,0,0,
              0,0,1,0,0,
              0,0,0,1,0,
              0,0,0,0,1,
              1,0,0,0,0), nrow = V, byrow = T)

Omega <- matrix(1, nrow = V, ncol = V)
Omega[1,3] <- Omega[3,1] <- 0

B.true <- matrix(rnorm(V^2), nrow = V) * B
Omega.true <- ifelse(Omega,runif(V^2, -1, 1),0) + diag(rep(4, V))
Omega.true[lower.tri(Omega.true,diag=T)] = 0
Omega.true = Omega.true + t(Omega.true)
Omega.true = Omega.true + diag(rowSums(abs(Omega.true)) + rchisq(dim(Omega.true)[1], df = 1))

sigma = solve(diag(rep(1, V)) - B.true) %*% Omega.true %*% t(solve(diag(rep(1, V)) - B.true))
sum(eigen(sigma)$values < 0)

Y <- t(mvrnorm(n = n, mu = rep(0, V), Sigma = sigma))
B.initial <- matrix(rnorm(V^2), nrow = V) * B
Omega.initial <- matrix(rnorm(V^2), nrow = V) %*% t(matrix(rnorm(V^2), nrow = V))
Omega.initial <- Omega.initial * Omega

sum((B.initial - B.true)^2)
sum((Omega.initial - Omega.true)^2)



sem_el_fitC(y_r = Y, b_r = B, omega_r = Omega,
            b_weights_r = B.true, d_r = rep(n, n),
            lambda_r = rep(5, V), gamma_r = rep(5, sum(Omega == 0)/2), v = V, tol = 1e-6, max_iter = 100)




