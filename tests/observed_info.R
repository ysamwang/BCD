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
sim.size <- 200
n.list <- c(100, 200,500, 1000, 5000)

for(type in c("gauss")){ 
  
  ci.record <- matrix(0, nrow = length(n.list), ncol = (sum(mod$B) + (sum(mod$Omega) - V)/2 + V) * 2)
  se.record <- matrix(0, nrow = length(n.list), ncol = (sum(mod$B) + (sum(mod$Omega) - V)/2 + V) * 2)
  entropy.record <- rep(0, length(n.list))
  for(n.index in 1:length(n.list)){
    b.mse.like <- b.mse.el <- o.mse.like <- o.mse.el <- b.bias.like <- b.bias.el <- o.bias.like <- o.bias.el <- rep(0, sim.size)
    ci.like <- ci.sw <- se.length.like <- se.length.sw <- rep(0, sum(mod$B) + (sum(mod$Omega) - V)/2 + V)
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
      
      
      
    I_obs <- getRicfInfo(S = out_ricf$SigmaHat, B = mod$B, Omega = mod$Omega ,
                       B.hat = out_ricf$BHat, Omega.hat = out_ricf$OmegaHat, type = "observed")
    sw <- getRicfSandwich(Y, mod$B, mod$Omega, out_ricf$BHat, out_ricf$OmegaHat)
    
      
      beta.like <- out_ricf$BHat[which(mod$B==1)]
      omega.like <- out_ricf$OmegaHat[lower.tri(mod$Omega.true, diag = T) & mod$Omega ==1]
      
      b.mse.like[i] <-  mean((beta.like - b.true)^2)
      b.bias.like[i] <- mean(beta.like - b.true)
      
      o.mse.like[i] <-  mean((omega.like - omega.true)^2)
      o.bias.like[i] <- mean(omega.like - omega.true)
      
      se.like <- sqrt(diag(solve(I_obs))/n)
      se.sw <- sqrt(diag(sw)/n)
      se.length.like <- se.length.like + se.like
      se.length.sw <- se.length.sw + se.sw
      ci.like <- ci.like +  ((abs(c(beta.like - b.true, omega.like - omega.true))  / se.like) < 1.959964)
      ci.sw <- ci.sw +  ((abs(c(beta.like - b.true, omega.like - omega.true))  / se.sw) < 1.959964)
    }
    ci.record[n.index, ] <- c(ci.like, ci.sw)/ sim.size
    se.record[n.index, ] <- c(se.length.like, se.length.sw)/ sim.size
  }
  saveRDS(ci.record, paste("ci_record_", type, ".RDS", sep = ""))
  saveRDS(se.record, paste("se_record_", type, ".RDS", sep = ""))
}
