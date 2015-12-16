library('BCD')
library('MASS')
library('mvtnorm')
library('batch')



set.seed(123)
V <- 5
n <- 1000
k <- 3
b <- .3
d <- .3
mod <- generateModel(v = V, n = n, k = k, b = b, d = d)
seed <- 1

Sigma <- solve(diag(rep(1, V)) - mod$B.true) %*% mod$Omega.true %*% t(solve(diag(rep(1, V)) - mod$B.true))

mult.factor <- solve(diag(rep(1,V)) - mod$B.true) %*% t(chol(mod$Omega.true))

###
parseCommandArgs()
set.seed(seed)


Omega.null <- mod$Omega
Omega.alt <- mod$Omega
Omega.alt[3,1] <- Omega.alt[1,3] <-1

    mu <- 0
    sigma <- 1 
    eps <- matrix(exp(rnorm(n * V, mean = mu, sd = sigma)),
                  nrow = V, ncol = n) - exp(mu + sigma^2/2)
    eps <- eps / sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2))
    Y <- mult.factor %*% eps
    Y <- Y - rowMeans(Y)

  
  ### RICF Procedure ###
  out_ricf_null <- ricf(mod$B, Omega.null, Y)
  out_ricf_alt <- ricf(mod$B, Omega.alt, Y)
  
  ricf.null <- -2 * (sum(dmvnorm(t(Y), sigma = out_ricf_null$SigmaHat, log = T))
                        - sum(dmvnorm(t(Y), sigma = out_ricf_alt$SigmaHat, log = T)))
  
  ### EL Procedure ###
  el_fit_null <- fitEL(Y, mod$B, Omega.null)
  el_fit_alt <- fitEL(Y, mod$B, Omega.alt)
  
  #       chol.sigma.el <- t(chol(Y %*% diag(el_fit_null$p) %*% t(Y)))
  #       chol.obs <- t(chol(Y %*% t(Y)))
  #       Y.null <- chol.sigma.el %*% solve(chol.obs, Y)
  #       for(k in 1:bartlett.boots){
  #         Y.b <- Y.null[, sample(n, size = n, replace = T)]
  #         el_fit_null.b <- fitEL(Y.b, mod$B, Omega.null)
  # 
  #         null.bartlett <- null.bartlett + -2 * sum(log(el_fit_null.b$p ))
  #         el_fit_alt.b <- fitEL(Y.b, mod$B, Omega.alt)
  # 
  #         alt.bartlett <- alt.bartlett + -2 * sum(log(el_fit_alt.b$p))
  #       }
  
  original.null <- -2 * sum(log(el_fit_null$p * n))
  #       corrected.null[i] <- -2 * sum(log(el_fit_null$p)) * 8 / 
  #         (null.bartlett / bartlett.boots)
  
  ### EL Procedure ###
  
  original.alt <- -2 * sum(log(el_fit_alt$p  * n))
  #       corrected.alt[i] <- -2 * sum(log(el_fit_alt$p )) * 7 / 
  #         (alt.bartlett / bartlett.boots)

result <- data.frame(original.null = original.null,
                     original.alt = original.alt,
                     ricf.null = ricf.null,
                     seed = seed)  
result_name = paste("~/BCD/results/res",seed, ".csv", sep = "")
write.csv(result, result_name, row.names = F)