set.seed(111)
library(MASS)
V <- 5
n <- 100
Omega <- matrix(0, nrow = V, ncol = V)
Omega[lower.tri(Omega)] <- rbinom(n = V*(V-1)/2, size = 1, prob = .3)
Omega.weights <- matrix(runif(V^2), nrow = V)
Omega.weights <- Omega.weights %*% t(Omega.weights)

Omega <- Omega + diag(rep(1,V)) + t(Omega)
Omega.weights <- Omega * Omega.weights

B <- matrix(0,nrow = V, ncol = V)
B[which(Omega==0)] <- rbinom(n = V^2-sum(Omega), size = 1, prob = .2)
B.weights <- matrix(rnorm(V^2), nrow = V) * B

sigma <- solve(diag(rep(1,V)) - B.weights) %*% Omega.weights %*% t(solve(diag(rep(1,V)) - B.weights))
ci_el <- rep(0,sum(B))
ci_ricf <- rep(0, sum(B))
init_b <- init_om <- el_b <- el_om <- ricf_b <- ricf_om <- 0


sim.size <- 100

for(i in 1:sim.size){
  
  ## Generate Data
  Y <- t(mvrnorm(n = n, mu = rep(0, V), Sigma = sigma))
  Y <- Y - rowMeans(Y)
  
  
  
  ## EL Routine 
  out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                        OmegaInit = NULL, sigConv = 0, maxIter = 1,
                        msgs = FALSE, omegaInitScale = .9)
  
  
  num_dual_vars = V + sum(Omega == 0)/2
  init_val = c(out_ricf_init$BHat[B==1])
  
  optim_out <- optim(par = init_val, fn = sem_el_fit_obj, gr = NULL,
                     y_r = Y, omega_r = Omega, b_r = B,
                     dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 60, meanEst = 0,
                     method = "BFGS", control = list(fnscale = -1))
  
  fitted_mod <- sem_el_fit_weights(optim_out$par, y_r = Y, omega_r = Omega, b_r = B,
                                   dual_r = rep(0, num_dual_vars), tol = 1e-14, max_iter = 100, meanEst = 0)
  
  B.hat <- matrix(0, nrow = V, ncol = V)
  B.hat[which(B==1)] <- optim_out$par
  
  Omega.hat <- (((diag(rep(1,V)) - B.hat) %*% Y)) %*% diag(c(1/fitted_mod$d)) %*% t((((diag(rep(1,V)) - B.hat) %*% Y)))
  
  sw <- empLike_sandwich(Y, B, Omega, B.hat, c(1/fitted_mod$d))
  se.el <- sqrt(diag(sw$sw.var)/n) 
  
  cover <- c(optim_out$par - B.weights[which(B==1)])/se.el
  
  ci_el <- ci_el + (abs(cover) < qnorm(.975))
  
  el_om <- el_om + mean(((Omega.hat - Omega.weights) * Omega)^2)
  el_b <- el_b + mean(((B.hat - B.weights) * B)^2)
  
  init_om <- init_om + mean(((out_ricf_init$OmegaHat - Omega.weights) * Omega)^2)
  init_b <- init_b + mean(((out_ricf_init$BHat - B.weights) * B)^2)
  
  
  
  ## RICF Procedure
  out_ricf <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                        OmegaInit = NULL, sigConv = 0, maxIter = 5000,
                        msgs = FALSE, omegaInitScale = .9)
  
  I <- ricf_info(sigma, Omega.weights, out_ricf$BHat, B)
  se.like <- sqrt(diag(solve(I))/n) 
  
  cover <- c(out_ricf$BHat[which(B==1)] - B.weights[which(B==1)])/se.like
  
  ci_ricf <- ci_ricf + (abs(cover) < qnorm(.975))
  
}

ci_el / sim.size
ci_ricf / sim.size

el_om / sim.size
el_b / sim.size

init_om / sim.size
init_b / sim.size
