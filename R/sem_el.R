sem_el <- function(b_par, y, omega, b)
{
  return(sem_el_fitC(y_r = Y, b_r = b, omega_r = omega,
              b_weights_r = b, dual_r = rep(0, V + sum(Omega == 0)/2),
              v = dim(b)[1], tol = 1e-6, max_iter = 50))
}
