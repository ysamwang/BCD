library('BCD')
library('MASS')
# setwd("~/BCD")
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
n.list <- c(100, 200,500)
boot.size <- 1000


type <- "lognormal"
record <- data.frame(b.mse.like = rep(0, length(n.list)),
                     b.bias.like = rep(0, length(n.list)),
                     o.mse.like = rep(0, length(n.list)),
                     o.bias.like = rep(0, length(n.list)))

ci.record <- matrix(0, nrow = length(n.list), ncol = (sum(mod$B) + (sum(mod$Omega) - V)/2 + V) * 2)
se.record <- matrix(0, nrow = length(n.list), ncol = (sum(mod$B) + (sum(mod$Omega) - V)/2 + V) * 2)
se.boot.temp <- matrix(0, nrow = boot.size, ncol = (sum(mod$B) + (sum(mod$Omega) - V)/2 + V))

entropy.record <- rep(0, length(n.list))
for(n.index in 1:length(n.list)){
  b.mse.like <- b.mse.el <- o.mse.like <- o.mse.el <- b.bias.like <- b.bias.el <- o.bias.like <- o.bias.el <- rep(0, sim.size)
  ci.like <- ci.boot <- se.length.like <- se.length.boot <- rep(0, sum(mod$B) + (sum(mod$Omega) - V)/2 + V)
  entropy <- 0
  n <- n.list[n.index]
  cat(paste("\n===",n, type,"\n")) 
  for(i in 1:sim.size){ 
    cat(paste(i,": ")) 
    
    mu <- 0
    sigma <- 1 
    eps <- matrix(exp(rnorm(n * V, mean = mu, sd = sigma)),
                  nrow = V, ncol = n) - exp(mu + sigma^2/2)
    eps <- eps / sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2))
    Y <- mult.factor %*% eps
    
    
    
    
    #### RICF Procedure ####
    out_ricf <- ricf(B = mod$B, Omega = mod$Omega, Y = Y, BInit = NULL,
                     OmegaInit = NULL, sigConv = 0, maxIter = 10000,
                     msgs = FALSE, omegaInitScale = .9, tol = 1e-6)
    
    
    
    I <- ricf_info(S = out_ricf$SigmaHat, BHat = out_ricf$BHat, OmegaHat = out_ricf$OmegaHat,
                   B = mod$B, Omega = mod$Omega)
    
    beta.like <- out_ricf$BHat[which(mod$B==1)]
    omega.like <- out_ricf$OmegaHat[lower.tri(mod$Omega.true, diag = T) & mod$Omega ==1]
    
    b.mse.like[i] <-  mean((beta.like - b.true)^2)
    b.bias.like[i] <- mean(beta.like - b.true)
    
    o.mse.like[i] <-  mean((omega.like - omega.true)^2)
    o.bias.like[i] <- mean(omega.like - omega.true)
    
    se.like <- sqrt(diag(solve(I))/n)
    se.length.like <- se.length.like + se.like 
    ci.like <- ci.like +  ((abs(c(beta.like - b.true, omega.like - omega.true))  / se.like) < 1.959964)
    
    
    for(i in 1:boot.size){
      Y.mod <- Y[,sample(n, replace = T)]
      out_ricf.boot <- ricf(B = mod$B, Omega = mod$Omega, Y = Y.mod, BInit = NULL,
                            OmegaInit = NULL, sigConv = 0, maxIter = 10000,
                            msgs = FALSE, omegaInitScale = .9, tol = 1e-6)
      
      se.boot.temp[i, ] <- c(out_ricf.boot$BHat[which(mod$B==1)],
                             out_ricf.boot$OmegaHat[lower.tri(mod$Omega.true, diag = T) & mod$Omega ==1])
    }
    se.boot.temp <- apply(se.boot.temp, 2, sort, decreasing = T)
    se.boot.high <- se.boot.temp[round(boot.size * .025), ]
    se.boot.low <- se.boot.temp[round(boot.size * .975), ]
    se.length.boot <- se.length.boot + (se.boot.high - se.boot.low) 
    
    ci.boot <- ci.boot +  
      (((se.boot.high - c(b.true, omega.true)) * (c(b.true, omega.true) - se.boot.low)) >0)
  } 
  record$b.mse.like[n.index] <- mean(b.mse.like)
  record$b.bias.like[n.index] <- mean(b.bias.like)
  record$o.mse.like[n.index] <- mean(o.mse.like)
  record$o.bias.like[n.index] <- mean(o.bias.like)
  ci.record[n.index, ] <- c(ci.like, ci.boot)/ sim.size
  se.record[n.index, ] <- c(se.length.like, se.length.boot)/ sim.size
}

saveRDS(ci.record, "bootstrap_record.RDS")
png("bootstrap_omega.png")
plot(ci.record[,31], ylim = c(.5, 1),
     main = "Bootstrapped CI's for Omega", type = "b", xaxt = "n",
     ylab = "Coverage Ratio", xlab = "Sample Size")
mtext("200 Simulation Runs; 1000 Bootstrapped Draws")
for(i in 32:46){
  lines(ci.record[, i], type = "b")
}
axis(1, at = c(1:3), labels = c(100,200,500))
abline(h = .95, lwd = 3, col = "red")
dev.off()