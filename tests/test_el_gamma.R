#### Testing Empirical Likelihood for Structural Equation Models ####
source("tests/simHelper.R")
library(microbenchmark)

set.seed(101)
V <- 10
n <- 200
k <- 5
b <- .2
d <- .2

sim.size <- 200
record <- data.frame(times = rep(0, sim.size),
                     times_naive = rep(0, sim.size),
                     times_euclid = rep(0, sim.size),
                     error_b = rep(0, sim.size),
                     error_b_naive = rep(0, sim.size),
                     error_b_euclid =rep(0, sim.size),
                     error_b_ricf =rep(0, sim.size),
                     error_omega = rep(0, sim.size),
                     error_omega_naive = rep(0, sim.size),
                     error_omega_euclid=rep(0, sim.size),
                     error_omega_ricf=rep(0, sim.size))

for(sim in 1:sim.size)
{
  cat(paste(sim, ". "))
  mod <- generateModel(v = V, n = n, k = k, b = b, d = d, errorDist = "pois")
  
  out_ricf <- ricf(B = mod$B, Omega = mod$Omega, Y = mod$Y, BInit = NULL,
                               OmegaInit = NULL, sigConv = 0, maxIter = 500,
                               msgs = FALSE, omegaInitScale = .9)
  
  
  num_dual_vars = V + sum(mod$Omega == 0)/2
  init_val = c(out_ricf$BHat[mod$B==1])
  
  # fitted_init <- sem_el_fit_weights(init_val,y_r = mod$Y, omega_r = mod$Omega, b_r = mod$B,
  #                    dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 100)
  
  mb <- microbenchmark(optim_out <- optim(par = init_val, fn = sem_el_fit_obj, gr = NULL,
        y_r = mod$Y, omega_r = mod$Omega, b_r = mod$B,
        dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 100,
        method = "BFGS", control = list(fnscale = -1))
        , times = 1, control = list(warmpup = 1))
  
  
  fitted_mod <- sem_el_fit_weights(optim_out$par,y_r = mod$Y, omega_r = mod$Omega, b_r = mod$B,
                                   dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 100)
  B <- matrix(0, nrow = V, ncol = V)
  B[mod$B==1] <- optim_out$par
  omega <- (diag(rep(1, V)) - B) %*% mod$Y %*% diag(c(1/fitted_mod$d)) %*% t(mod$Y) %*% (diag(rep(1, V)) - B)
  
  
  mb_euclid <- microbenchmark(optim_out_euclid <- optim(par = init_val, fn = sem_el_euclid_fit_obj, gr = NULL,
                                          y_r = mod$Y, omega_r = mod$Omega, b_r = mod$B,
                                          method = "BFGS", control = list(fnscale = -1))
                       , times = 1, control = list(warmpup = 1))
  
  fitted_mod_euclid <- sem_el_euclid_fit_weights(optim_out_euclid$par,y_r = mod$Y, omega_r = mod$Omega, b_r = mod$B)
  B <- matrix(0, nrow = V, ncol = V)
  B[mod$B==1] <- optim_out_euclid$par
  omega_euclid <- (diag(rep(1, V)) - B) %*% mod$Y %*% diag(c(fitted_mod_euclid$p_star)) %*% t(mod$Y) %*% t((diag(rep(1, V)) - B))
  
  omeg_weights <- c()
  omeg_truth <- c()
  omeg_euclid <- c()
  omeg <- c()
  omeg_ricf <- c()
  for(i in 1:V){
    for(j in i:V){
      if(mod$Omega[i,j] == 1)
      {
        omeg_weights <- c(omeg_weights, out_ricf$OmegaHat[i,j])
        omeg_truth <- c(omeg_truth, mod$Omega.true[i,j])
        omeg_euclid <- c(omeg_euclid, omega_euclid[i,j])
        omeg <- c(omeg, omega[i,j])
        omeg_ricf <- c(omeg_ricf, out_ricf$OmegaHat[i,j])
      }
    }
  }
  
  init_val_naive = c(out_ricf$BHat[mod$B==1], omeg_weights)
  
  num_dual_vars_naive = V + V * (V + 1) / 2
  mb_naive <- microbenchmark(optim_out_naive <- optim(par = init_val_naive, fn = sem_el_naive_fit_obj, gr = NULL,
               y_r = mod$Y, omega_r = mod$Omega, b_r = mod$B,
               dual_r = rep(0, num_dual_vars_naive), tol = 1e-6, max_iter = 100,
               method = "BFGS", control = list(fnscale = -1)),
               times = 1, control = list(warmpup = 1))
  
  fitted_mod_naive <- sem_el_naive_fit_weights(optim_out_naive$par,y_r = mod$Y, omega_r = mod$Omega, b_r = mod$B,
                                   dual_r = rep(0, num_dual_vars_naive), tol = 1e-6, max_iter = 100)

  record$times[sim] <- mb$time * (optim_out$convergence == 0) * (sum(1/fitted_mod$d) == 1)
  record$times_naive[sim] <- mb_naive$time  * (optim_out_naive$convergence == 0) * (sum(1/fitted_mod_naive$d) == 1)
  record$times_euclid[sim] <- mb_euclid$time * (optim_out_euclid$convergence == 0) * (sum(fitted_mod_euclid$p_star) == 1)
  
  record$error_b_ricf[sim] <- sum((c(out_ricf$BHat[which(c(mod$B)==1)]) - c(mod$B.true[which(c(mod$B)==1)]))^2)
  record$error_b[sim] <- sum((optim_out$par - c(mod$B.true[which(c(mod$B)==1)]))^2)
  record$error_b_naive[sim] <- sum((optim_out_naive$par[1:sum(mod$B)] - c(mod$B.true[which(c(mod$B)==1)]))^2)
  record$error_b_euclid[sim] <-sum((optim_out_euclid$par - c(mod$B.true[which(c(mod$B)==1)]))^2)
  
  record$error_omega_ricf[sim] <- sum((omeg_ricf - omeg_truth)^2)
  record$error_omega[sim] <- sum((omeg_ricf - omeg_truth)^2)
  record$error_omega_naive[sim] <- sum((optim_out_naive$par[(sum(mod$B)+1):(sum(mod$B)+length(omeg_weights))] - omeg_truth)^2)
  record$error_omega_euclid[sim] <- sum((omeg_euclid - omeg_truth)^2)  
}
