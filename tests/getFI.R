source("tests/simHelper.R")

set.seed(123)
V <- 8
n <- 2
k <- 3
b <- .3
d <- .3
dist <- "gauss"
mod <- generateModel(v = V, n = n, k = k, b = b, d = d, errorDist = dist)

Sigma <- solve(diag(rep(1, V)) - mod$B.true) %*% mod$Omega.true %*% t(solve(diag(rep(1, V)) - mod$B.true))

mult.factor <- solve(diag(rep(1,V)) - mod$B.true) %*% t(chol(mod$Omega.true))

omega.true <- mod$Omega.true[lower.tri(mod$Omega.true, diag = T)][which(mod$Omega[lower.tri(mod$Omega, diag = T)]==1)]
sim.size <- 1000
n.list <- c(100, 200,500, 1000, 2000)

for(type in c("gauss", "lognormal")){ 
  record <- data.frame(mse.like = rep(0, length(n.list)),
                       bias.like = rep(0, length(n.list)),
                       mse.el = rep(0, length(n.list)),
                       bias.el = rep(0, length(n.list)))
  
  ci.record <- matrix(0, nrow = length(n.list), ncol = sum(mod$B) * 2)
  se.record <- matrix(0, nrow = length(n.list), ncol = sum(mod$B) * 2)
  entropy.record <- rep(0, length(n.list))
  
  for(n.index in 1:length(n.list)){
    mse.like <- mse.el <- mse.init <- bias.like <- bias.init <- bias.el <- rep(0, sim.size)
    ci.like <- ci.el <- se.length.like <- se.length.el <- rep(0, sum(mod$B))
    entropy <- 0
    for(i in 1:sim.size){
      cat(paste(i,": "))
      n <- n.list[n.index]
      
      if(type == "gauss"){
        Y <- t(MASS::mvrnorm(n = n, Sigma = Sigma, mu = rep(0, V)))
        Y <- Y - rowMeans(Y)
      } else {
        mu <- 0
        sigma <- 1 
        eps <- matrix(exp(rnorm(n * V, mean = mu, sd = sigma)),
                      nrow = V, ncol = n) - exp(mu + sigma^2/2)
        eps <- eps / sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2))
        Y <- mult.factor %*% eps
      }
      
      
      
      
      #### RICF Procedure ####
      out_ricf <- ricf(B = mod$B, Omega = mod$Omega, Y = Y, BInit = NULL,
                                                 OmegaInit = NULL, sigConv = 0, maxIter = 10000,
                                                 msgs = FALSE, omegaInitScale = .9, tol = 1e-6)
      
      
      
      I <- ricf_info(S = out_ricf$SigmaHat, BHat = out_ricf$BHat, OmegaHat = out_ricf$OmegaHat,
                     B = mod$B, Omega = mod$Omega)
      
      beta.like <- out_ricf$BHat[which(mod$B==1)]
      mse.like[i] <-  mean((beta.like - mod$B.true[which(mod$B==1)])^2)
      bias.like[i] <- mean(beta.like - mod$B.true[which(mod$B==1)])
      
#       omega.like <- out_ricf$OmegaHat[lower.tri(out_ricf$OmegaHat, diag = T)][which(mod$Omega[lower.tri(mod$Omega, diag = T)]==1)]
#       mse.like[i] <-  mean((omega.like - omega.true)^2)
#       bias.like[i] <- mean(omega.like - omega.true)

            
#       se.like <- sqrt(diag(solve(I))/n)[-c(1:sum(mod$B))]
#       se.length.like <- se.length.like + se.like 
#       ci.like <- ci.like +  (abs(omega.like - omega.true) / se.like < qnorm(.975))
      
      se.like <- sqrt(diag(solve(I))/n)[c(1:sum(mod$B))]
      se.length.like <- se.length.like + se.like 
      ci.like <- ci.like +  (abs(beta.like - mod$B.true[which(mod$B==1)]) / se.like < qnorm(.975))
    
      ##### EL Procedure #####
      
      out_ricf_init <- ricf(B = mod$B, Omega = mod$Omega, Y = Y, BInit = NULL,
                            OmegaInit = NULL, sigConv = 0, maxIter = 1,
                            msgs = FALSE, omegaInitScale = .9)
      
      
      num_dual_vars = V + sum(mod$Omega == 0)/2
      init_val = c(out_ricf_init$BHat[mod$B==1])
      mse.init[i] <- mean((init_val - mod$B.true[which(mod$B==1)])^2)
      bias.init[i] <- mean(init_val - mod$B.true[which(mod$B==1)])
      
      optim_out <- optim(par = init_val, fn = sem_el_fit_obj, gr = NULL,
                         y_r = Y, omega_r = mod$Omega, b_r = mod$B,
                         dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 60,
                         meanEst = 0,
                         method = "BFGS", control = list(fnscale = -1))
      
      fitted_mod <- sem_el_fit_weights(optim_out$par,y_r = Y, omega_r = mod$Omega, b_r = mod$B,
                                       dual_r = rep(0, num_dual_vars), tol = 1e-6, max_iter = 100, meanEst = 0)
      B.hat <- matrix(0, nrow = V, ncol = V)
      B.hat[mod$B==1] <- optim_out$par
      
      sw <- empLike_sandwich(Y, mod$B, mod$Omega, B.hat, c(1/fitted_mod$d))
      
      beta.el <- B.hat[which(mod$B==1)]
      mse.el[i] <- mean((optim_out$par - mod$B.true[which(mod$B==1)])^2)
      bias.el[i] <- mean(optim_out$par - mod$B.true[which(mod$B==1)])
      
      se.el <- sqrt(diag(sw$sw.var)/n)
      se.length.el <- se.length.el + se.el
      ci.el <- ci.el +  (abs(beta.el - mod$B.true[which(mod$B==1)]) / se.el < qnorm(.975))
      entropy <- entropy + -sum(c(1/fitted_mod$d) * log(c(1/fitted_mod$d)))
    } 
    record$mse.like[n.index] <- mean(mse.like)
    record$bias.like[n.index] <- mean(bias.like)
    record$mse.el[n.index] <- mean(mse.el)
    record$bias.el[n.index] <- mean(bias.el)
    ci.record[n.index, ] <- c(ci.like, ci.el)/ sim.size
    se.record[n.index, ] <- c(se.length.like, se.length.el)/ sim.size
    entropy.record[n.index] <- entropy / sim.size
  }
  saveRDS(record, paste("record_", type, sep = ""))
  saveRDS(ci.record, paste("ci_record_", type, sep = ""))
  saveRDS(se.record, paste("se_record_", type, sep = ""))
  saveRDS(entropy.record, paste("entropy_record_", type, sep = ""))
}

