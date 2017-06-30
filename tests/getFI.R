library('BCD')
library('MASS')

set.seed(123)
V <- 8
n <- 1000
k <- 3
b <- .3
d <- .3
dist <- "gauss"
mod <- generateModel(v = V, n = n, k = k, b = b, d = d, errorDist = dist)

Sigma <- solve(diag(rep(1, V)) - mod$B.true) %*% mod$Omega.true %*% t(solve(diag(rep(1, V)) - mod$B.true))

mult.factor <- solve(diag(rep(1,V)) - mod$B.true) %*% t(chol(mod$Omega.true))

omega.true <- mod$Omega.true[lower.tri(mod$Omega.true, diag = T)][which(mod$Omega[lower.tri(mod$Omega, diag = T)]==1)]
b.true <- mod$B.true[which(mod$B==1)]
sim.size <- 500
n.list <- c(100, 200,500, 1000, 5000)

for(type in c("gauss", "lognormal")){ 
  record <- data.frame(b.mse.like = rep(0, length(n.list)),
                       b.bias.like = rep(0, length(n.list)),
                       b.mse.el = rep(0, length(n.list)),
                       b.bias.el = rep(0, length(n.list)),
                       o.mse.like = rep(0, length(n.list)),
                       o.bias.like = rep(0, length(n.list)),
                       o.mse.el = rep(0, length(n.list)),
                       o.bias.el = rep(0, length(n.list)))
  
  ci.record <- matrix(0, nrow = length(n.list), ncol = (sum(mod$B) + (sum(mod$Omega) - V)/2 + V) * 3)
  se.record <- matrix(0, nrow = length(n.list), ncol = (sum(mod$B) + (sum(mod$Omega) - V)/2 + V) * 3)
  entropy.record <- rep(0, length(n.list))
  for(n.index in 1:length(n.list)){
    b.mse.like <- b.mse.el <- o.mse.like <- o.mse.el <- b.bias.like <- b.bias.el <- o.bias.like <- o.bias.el <- rep(0, sim.size)
    ci.like <- ci.el <- ci.like_obs <- se.length.like <- se.length.el <- se.length.like_obs <- rep(0, sum(mod$B) + (sum(mod$Omega) - V)/2 + V)
    entropy <- 0
    n <- n.list[n.index]
    cat(paste("\n===",n, type,"\n"))
    for(i in 1:sim.size){
      cat(paste(i,": "))

      
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
      
      
      
      I <- getRicfInfo(S = out_ricf$SigmaHat, B = mod$B, Omega = mod$Omega ,B.hat = out_ricf$BHat, Omega.hat = out_ricf$OmegaHat,type = "expected")
      
      beta.like <- out_ricf$BHat[which(mod$B==1)]
      omega.like <- out_ricf$OmegaHat[lower.tri(mod$Omega.true, diag = T) & mod$Omega ==1]
      
      b.mse.like[i] <-  mean((beta.like - b.true)^2)
      b.bias.like[i] <- mean(beta.like - b.true)
      
      o.mse.like[i] <-  mean((omega.like - omega.true)^2)
      o.bias.like[i] <- mean(omega.like - omega.true)
      
      se.like <- sqrt(diag(solve(I))/n)
      se.length.like <- se.length.like + se.like 
      ci.like <- ci.like +  ((abs(c(beta.like - b.true, omega.like - omega.true))  / se.like) < 1.959964)
  
      ci.like <- ci.like +  ((abs(c(beta.like - b.true, omega.like - omega.true))  / se.like) < 1.959964)  
      ##### EL Procedure #####
      
      out_ricf_init <- ricf(B = mod$B, Omega = mod$Omega, Y = Y, BInit = NULL,
                            OmegaInit = NULL, sigConv = 0, maxIter = 1,
                            msgs = FALSE, omegaInitScale = .8)
      
      
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
      omega.el <- Omega.hat[lower.tri(Omega.hat, diag = TRUE) & mod$Omega ==1]
      
      sw <- empLike_sandwich(Y, mod$B, mod$Omega, B.hat, Omega.hat, c(1/fitted_mod$d))
      
      beta.el <- optim_out$par
      b.mse.el[i] <- mean((beta.el - b.true)^2)
      b.bias.el[i] <- mean(beta.el - b.true)
      
      o.mse.el[i] <- mean((omega.el - omega.true)^2)
      o.bias.el[i] <- mean(omega.el - omega.true)
      
      se.el <- sqrt(diag(sw$sw.var)/n)
      se.length.el <- se.length.el + se.el
      ci.el <- ci.el +  ((abs(c(beta.el - b.true,  omega.el - omega.true))  / se.el) < 1.959964)
      entropy <- entropy + -sum(c(1/fitted_mod$d) * log(c(1/fitted_mod$d)))
    } 
    record$b.mse.like[n.index] <- mean(b.mse.like)
    record$b.bias.like[n.index] <- mean(b.bias.like)
    record$b.mse.el[n.index] <- mean(b.mse.el)
    record$b.bias.el[n.index] <- mean(b.bias.el)
    
    record$o.mse.like[n.index] <- mean(o.mse.like)
    record$o.bias.like[n.index] <- mean(o.bias.like)
    record$o.mse.el[n.index] <- mean(o.mse.el)
    record$o.bias.el[n.index] <- mean(o.bias.el)
    
    
    ci.record[n.index, ] <- c(ci.like, ci.el)/ sim.size
    se.record[n.index, ] <- c(se.length.like, se.length.el)/ sim.size
    entropy.record[n.index] <- entropy / sim.size
  }
  saveRDS(record, paste("record_", type, ".RDS", sep = ""))
  saveRDS(ci.record, paste("ci_record_", type, ".RDS", sep = ""))
  saveRDS(se.record, paste("se_record_", type, ".RDS", sep = ""))
  saveRDS(entropy.record, paste("entropy_record_", type,".RDS", sep = ""))
}
