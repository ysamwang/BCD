#### BCD Test File ####

# True Values acyclic B
library(MASS)
library(microbenchmark)
set.seed(100)
V <- 20
n <- 750

B <- matrix(0, nrow = V, ncol = V)
for(i in 2:V)
{
  B[i-1, i] = 1
}
B[V, 1] = 1

Omega <- matrix(0, nrow = V, ncol = V)
Omega[upper.tri(Omega, diag = F)] <- rbinom(size = 1, prob = .5, n = V * (V-1)/ 2)
Omega[c(B)== 1] <- 0
Omega = Omega + t(Omega) + diag(rep(1, V))

B.true <- matrix(rnorm(V^2, mean = 0, sd = 3), nrow = V) * B
Omega.true <- ifelse(Omega,runif(V^2, -1, 1),0) + diag(rep(4, V))
Omega.true[lower.tri(Omega.true,diag=T)] = 0
Omega.true = Omega.true + t(Omega.true)
Omega.true = Omega.true + diag(rowSums(abs(Omega.true)) + rchisq(dim(Omega.true)[1], df = 1))

sigma = solve(diag(rep(1, V)) - B.true) %*% Omega.true %*% t(solve(diag(rep(1, V)) - B.true))
sum(eigen(sigma)$values < 0)

Y <- t(mvrnorm(n = n, mu = rep(0, V), Sigma = sigma))
Y <- Y - rowMeans(Y)

#####

microbenchmark(sem_el_fit_obj(b_weights_r = c(B.true[B == 1]), y_r = Y, omega_r = Omega,
               b_r = B, dual_r = rep(0, V + sum(Omega == 0)/2), tol = 1e-6, max_iter = 100), times = 5)
omeg_weights <- c()
for(j in 1:V){
  for(i in j:V){
    if(Omega[i,j] ==1){
      omeg_weights <- c(omeg_weights, Omega.true[i,j])
    }
  }
}
microbenchmark(sem_el_naive_fit_obj(weights_r = c(B.true[B == 1], omeg_weights), y_r = Y, omega_r = Omega, 
                              b_r = B, dual_r = rep(0, V + V*(V+1)/2), tol = 1e-6, max_iter = 100), times = 5)
out <- sem_el_fit_weights(b_weights_r = c(B.true[B == 1]), y_r = Y, omega_r = Omega,
                      b_r = B, dual_r = rep(0, V + sum(Omega == 0)/2), tol = 1e-6, max_iter = 100)

out.naive <- sem_el_naive_fit_weights(weights_r = c(B.true[B == 1], omeg_weights), y_r = Y, omega_r = Omega, 
                     b_r = B, dual_r = rep(0, V + V*(V+1)/2), tol = 1e-6, max_iter = 100)
#####
st <- proc.time()
optim(c(B.true[B == 1]), sem_el_fit_obj, gr = NULL, y_r = Y, omega_r = Omega,
      b_r = B, dual_r = rep(0, V + sum(Omega == 0)/2), tol = 1e-6, max_iter = 100, method = "Nelder-Mead")
end <- proc.time()

st.naive <-proc.time()
optim(c(B.true[B == 1], omeg_weights), sem_el_fit_obj, gr = NULL, y_r = Y, omega_r = Omega, 
      b_r = B, dual_r = rep(0, V + V*(V+1)/2), tol = 1e-6, max_iter = 100, method = "Nelder-Mead")
end.naive <- proc.time()
