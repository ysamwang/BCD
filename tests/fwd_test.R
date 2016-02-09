set.seed(1000)

# Simulation Parameters
V <- 5
Omega <- diag(rep(1, V))
sim.size <- 100
dist <- "lognormal"
B.true.struct <- matrix(c(0, 0, 0, 0, 0,
                          1, 0, 0, 0, 0,
                          1, 0, 0, 0, 0,
                          0, 1, 0, 0, 0,
                          0, 1, 1, 0, 0), nrow = V, ncol = V, byrow = T)

n <- 500




### Pre-Processsing
B.true <- round(rnorm(V^2),2) * B.true.struct
Omega.true <- diag(rep(1, V))


temp <- solve(diag(rep(1,V)) - B.true)
sig.ln <- log(1 + (exp(1) - 1) * Omega.true)


sigma <- temp %*% Omega.true %*% t(temp)

errs <- matrix(0, nrow = 1, ncol = 2)
colnames(errs) <- c("FP", "FN")


## Actual Simulation
for(sim in 1:sim.size){
  cat("Sim: ")
  cat(sim)
  cat("\n")
  #### Generate Data ####
  if(dist == "t"){
    t.df <- 30
    errs <- rbind(rt(n, df = 5),
                  rt(n, df = 10),
                  rt(n, df = 15))
    Y <- temp %*% errs
    
  } else if(dist == "lognormal"){
    sig.ln <- log(1 + (exp(1) - 1) * Omega.true)
    z <- t(MASS::mvrnorm(n = n, mu = rep(0, V), Sigma = sig.ln))
    Y <- temp %*% ((exp(z) - exp(1/2))/ sqrt(exp(1) * (exp(1) - 1)))
  } else if(dist == "gauss") {
    Y <- temp %*% (t(MASS::mvrnorm(n = n, mu = rep(0, V), Sigma = Omega.true)))
  } else if(dist == "laplace"){
#     library(distr)
#     D <- DExp(rate = sqrt(2))
#     errs <-  t(chol(Omega.true)) %*% matrix(r(D)(n * V), nrow = V)
    Y <- temp %*% errs
  } else if(dist == "gamma"){
    errs <-  t(chol(Omega.true)) %*% matrix(rgamma(n * V, 1, 1) - 1, nrow = V)
    Y <- temp %*% errs
  } else if(dist == "unif"){
    errs <-  t(chol(Omega.true)) %*% matrix(sample((-2):(2), size = n *V, replace = T)/ sqrt(2), nrow = V)
    Y <- temp %*% errs
  }
  
output <- forwardBackwardDiscovery(Y, moments = 4, tol = 30, type = "profile")
  errs[sim, 1] <- sum(output == 1 & B.true.struct == 0)
  errs[sim, 2] <- sum(output == 0 & B.true.struct == 1)
}

