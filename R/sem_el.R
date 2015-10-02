sem_el <- function(b_par, y, omega, b, )
{
  return(sem_el_fitC(y_r = Y, omega_r = omega,
              b_weights_r = B.true, d_r = rep(n, n),
              dual_r = rep(0, V + sum(Omega == 0)/2), v = V, tol = 1e-6, max_iter = 50))
}