set.seed(1001)

V <- 3
n <- 300

B <- matrix(0, nrow = V, ncol = V)
Omega <- diag(rep(1,V))

sigma <- matrix(c(1, .5, .5,
                  .5, 1, .5,
                  .5, .5, 1), nrow = V, ncol = V)

Y <- t(MASS::mvrnorm(n = n, Sigma = sigma, mu = rep(0, V)))
Y  <- Y

#############
B <- matrix(0, nrow = V, ncol = V)
B[1,2] <- B[2,3]<- 1


out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                      OmegaInit = NULL, sigConv = 0, maxIter = 1,
                      msgs = FALSE, omegaInitScale = .9)


num_dual_vars = V + sum(Omega == 0)/2
init_val = c(rep(2,V), out_ricf_init$BHat[B==1])

optim_out1 <- optim(par = init_val, fn = sem_el_fit_obj, gr = NULL,
                   y_r = Y, omega_r = Omega, b_r = B,
                   dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 60,
                   method = "BFGS", control = list(fnscale = -1))

fitted_mod1 <- sem_el_fit_weights(optim_out1$par, y_r = Y, omega_r = Omega, b_r = B,
                                 dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 100)
optim_out1$par


#############
B <- matrix(0, nrow = V, ncol = V)
B[2,1] <- B[3,2] <-1


out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                      OmegaInit = NULL, sigConv = 0, maxIter = 1,
                      msgs = FALSE, omegaInitScale = .9)


num_dual_vars = V + sum(Omega == 0)/2
init_val = c(out_ricf_init$BHat[B==1])

optim_out2 <- optim(par = init_val, fn = sem_el_fit_obj, gr = NULL,
                    y_r = Y, omega_r = Omega, b_r = B,
                    dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 60,
                    method = "BFGS", control = list(fnscale = -1))

fitted_mod2 <- sem_el_fit_weights(optim_out2$par,y_r = Y, omega_r = Omega, b_r = B,
                                  dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 100)




#############
set.seed(888)
V <- 3
B <- matrix(0, nrow = V, ncol = V)
Omega <- diag(rep(1,V))

sigma <- matrix(c(1, .5, .5,
                  .5, 1, .5,
                  .5, .5, 1), nrow = V, ncol = V)

Y <- t(MASS::mvrnorm(n = n, Sigma = sigma, mu = rep(0, V)))
Y  <- Y - rowMeans(Y)

B <- matrix(0, nrow = V, ncol = V)
B[2,1] <- B[3,2] <-1


out_ricf_init <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                      OmegaInit = NULL, sigConv = 0, maxIter = 1,
                      msgs = FALSE, omegaInitScale = .9)


num_dual_vars = V + sum(Omega == 0)/2
init_val = c(out_ricf_init$BHat[B==1])

optim_outSEM <- optim(par = init_val, fn = sem_el_fit_obj, gr = NULL,
                    y_r = Y, omega_r = Omega, b_r = B,
                    dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 60,
                    method = "BFGS", control = list(fnscale = -1))

fitted_modSEM <- sem_el_fit_weights(optim_outSEM$par,y_r = Y, omega_r = Omega, b_r = B,
                                  dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 100)


##
V <- 2
Omega <- diag(rep(1,V))
B <- matrix(0, nrow = V, ncol = V)
B[2,1] <-1
Ymod <- Y[1:2,]


out_ricf_init <- ricf(B = B, Omega = Omega, Y = Ymod, BInit = NULL,
                      OmegaInit = NULL, sigConv = 0, maxIter = 1,
                      msgs = FALSE, omegaInitScale = .9)


num_dual_vars = V + sum(Omega == 0)/2
init_val = c(optim_outSEM$par[1])

optim_outReg21 <- optim(par = init_val, fn = sem_el_fit_obj, gr = NULL,
                     y_r = Ymod, omega_r = Omega, b_r = B,
                     dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 60,
                     method = "BFGS", control = list(fnscale = -1))

fitted_modReg21 <- sem_el_fit_weights(optim_outReg21$par,y_r = Ymod, omega_r = Omega, b_r = B,
                                    dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 100)


##
V <- 2
Omega <- diag(rep(1,V))
B <- matrix(0, nrow = V, ncol = V)
B[2,1] <-1
Ymod <- Y[2:3,]


out_ricf_init <- ricf(B = B, Omega = Omega, Y = Ymod, BInit = NULL,
                      OmegaInit = NULL, sigConv = 0, maxIter = 1,
                      msgs = FALSE, omegaInitScale = .9)


num_dual_vars = V + sum(Omega == 0)/2
init_val = c(optim_outSEM$par[2])

optim_outReg32 <- optim(par = init_val, fn = sem_el_fit_obj, gr = NULL,
                        y_r = Ymod, omega_r = Omega, b_r = B,
                        dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 60,
                        method = "BFGS", control = list(fnscale = -1))

fitted_modReg32 <- sem_el_fit_weights(optim_outReg32$par,y_r = Ymod, omega_r = Omega, b_r = B,
                                    dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 100)

