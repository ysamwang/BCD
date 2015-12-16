
.fit.single <- function(Y, B, Omega, B.hat, cutoff, row.ind, col.ind) {
  tol = 1e-6
  maxInnerIter = 100
  method = "BFGS"
  meanEst = 0  
  
  B.mod <- B 
  B.mod[row.ind, col.ind] <- 0
  
  init_val <- which(B.hat[B.mod == 1])
  
  
  optim_out <- optim(par = init_val, fn = sem_el_fit_obj_one_fixed, gr = NULL,
                     y_r = Y, omega_r = Omega, b_r = B.mod,
                     dual_r = rep(0, V + sum(Omega == 0)/2),
                     tol = tol, max_iter = maxInnerIter,
                     meanEst = meanEst, fixed = B.hat[row.ind, col.ind], row_ind = row.ind, col_ind = col.ind,
                     method = method, control = list(fnscale = -1))
  
  fitted_mod <- sem_el_fit_weights(optim_out$par,y_r = Y, omega_r = Omega, b_r = B,
                                   dual_r = rep(0, num_dual_vars), tol = tol, max_iter = maxInnerIter, meanEst = meanEst)
  
  p <- c(1/fitted_mod$d)
   
  # return LRT
  return(-2 * sum(log(p * n)))
}