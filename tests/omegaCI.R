source("tests/simHelper.R")

set.seed(123)
V <- 8
n <- 2000
k <- 3
b <- .3
d <- .3
dist <- "lognormal"

mod <- generateModel(v = V, n = n, k = k, b = b, d = d, errorDist = dist)
mult.factor <- solve(diag(rep(1,V)) - mod$B.true) %*% t(chol(mod$Omega.true))

omega.true <- mod$Omega.true[lower.tri(mod$Omega.true, diag = T)][which(mod$Omega[lower.tri(mod$Omega, diag = T)]==1)]
b.true <- mod$B.true[which(mod$B==1)]
sim.size <- 500

o.mse.like <- o.mse.el <- o.mse.oracle <- o.bias.like <- o.bias.el <- o.bias.oracle <-
  o.like <- o.el <- o.oracle <- matrix(0, nrow = sim.size, ncol = (sum(mod$Omega) - V) / 2 + V)

ci.like <- ci.el <- se.length.like <- se.length.el <- rep(0, sum(mod$B) + (sum(mod$Omega) - V)/2 + V)
      
for(i in 1:sim.size){
  cat(paste(i, ":"))
      mu <- 0
      sigma <- 1 
      eps <- matrix(exp(rnorm(n * V, mean = mu, sd = sigma)), nrow = V, ncol = n)
      eps <- eps - exp(mu + sigma^2/2)
      eps <- eps / sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2))
      eps <- t(chol(mod$Omega.true)) %*% eps
      Y <- solve(diag(rep(1,V)) - mod$B.true, eps)
    
    
    #### Baseline Measures
    o.oracle[i,] <- (eps %*% t(eps) / n)[lower.tri(mod$Omega, diag = T) & mod$Omega == 1]
    o.mse.oracle[i, ] <- (o.oracle[i,] - omega.true)^2
    o.bias.oracle[i, ] <- (o.oracle[i,] - omega.true)
    
    
    #### RICF Procedure ####
    out_ricf <- ricf(B = mod$B, Omega = mod$Omega, Y = Y, BInit = NULL,
                     OmegaInit = NULL, sigConv = 0, maxIter = 10000,
                     msgs = FALSE, omegaInitScale = .9, tol = 1e-6)
    
    
    
    I <- ricf_info(S = out_ricf$SigmaHat, BHat = out_ricf$BHat, OmegaHat = out_ricf$OmegaHat,
                   B = mod$B, Omega = mod$Omega)
    
    beta.like <- out_ricf$BHat[which(mod$B==1)]
    o.like[i,] <- out_ricf$OmegaHat[lower.tri(mod$Omega.true, diag = T) & mod$Omega ==1]
    
    o.mse.like[i,] <-  (o.like[i,] - omega.true)^2
    o.bias.like[i,] <- (o.like[i,] - omega.true)
    
    se.like <- sqrt(diag(solve(I))/n)
    se.length.like <- se.length.like + se.like 
    ci.like <- ci.like +  ((abs(c(beta.like - b.true, o.like[i,] - omega.true))  / se.like) < 1.959964)
    
    ##### EL Procedure #####
    
    out_ricf_init <- ricf(B = mod$B, Omega = mod$Omega, Y = Y, BInit = NULL,
                          OmegaInit = NULL, sigConv = 0, maxIter = 1,
                          msgs = FALSE, omegaInitScale = .9)
    
    
    num_dual_vars = V + sum(mod$Omega == 0)/2
    init_val = c(out_ricf_init$BHat[mod$B==1])
    
    optim_out <- optim(par = init_val, fn = sem_el_fit_obj, gr = NULL,
                       y_r = Y, omega_r = mod$Omega, b_r = mod$B,
                       dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 60,
                       meanEst = 0,
                       method = "BFGS", control = list(fnscale = -1))
    
    fitted_mod <- sem_el_fit_weights(optim_out$par,y_r = Y, omega_r = mod$Omega, b_r = mod$B,
                                     dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 60, meanEst = 0)
    
    B.hat <- matrix(0, nrow = V, ncol = V)
    B.hat[mod$B==1] <- optim_out$par
    
    Omega.hat <- (diag(rep(1,V)) - B.hat) %*% Y %*% diag(c(1/fitted_mod$d)) %*% t(Y) %*% t((diag(rep(1,V)) - B.hat))
    o.el[i,] <- Omega.hat[lower.tri(Omega.hat, diag = TRUE) & mod$Omega ==1]
    
    sw <- empLike_sandwich(Y, mod$B, mod$Omega, B.hat, Omega.hat, c(1/fitted_mod$d))  
    
    beta.el <- optim_out$par
    
    o.mse.el[i,] <- (o.el[i,] - omega.true)^2
    o.bias.el[i,] <- (o.el[i,] - omega.true)
    
    se.el <- sqrt(diag(sw$sw.var)/n)
    se.length.el <- se.length.el + se.el
    ci.el <- ci.el +  ((abs(c(beta.el - b.true,  o.el[i,] - omega.true))  / se.el) < 1.959964) 
}

saveRDS(list(o.oracle = o.oracle,
               o.like = o.like,
               o.el = o.el,
               o.mse.like = o.mse.like,
               o.mse.el = o.mse.el,
               o.mse.oracle = o.mse.oracle,
               o.bias.like = o.bias.like,
               o.bias.el = o.bias.el,
               o.bias.oracle = o.bias.oracle,
               se.length.like = se.length.like,
               se.length.el = se.length.el,
               ci.like = ci.like,
               ci.el = ci.el), "omegaCI.RDS")
