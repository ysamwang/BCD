set.seed(300)

V <- 5
B.true.struct <- matrix(c(0, 0, 0, 0, 0,
              1, 0, 0, 0, 0,
              1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 1, 1, 0, 0), nrow = V, ncol = V, byrow = T)
Omega <- diag(rep(1, V))



n <- 500


B.true <- round(rnorm(V^2),2) * B.true.struct
Omega.true <- diag(rep(1, V))


temp <- solve(diag(rep(1,V)) - B.true)
sig.ln <- log(1 + (exp(1) - 1) * Omega.true)


sigma <- temp %*% Omega.true %*% t(temp)
dist <- "lognormal"


sim.size <- 1
edit.dist <- rep(0, sim.size)
for(ind in 1:sim.size){
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
#   library(distr)
#   D <- DExp(rate = sqrt(2))
#   errs <-  t(chol(Omega.true)) %*% matrix(r(D)(n * V), nrow = V)
#   Y <- temp %*% errs
} else if(dist == "gamma"){
  errs <-  t(chol(Omega.true)) %*% matrix(rgamma(n * V, 1, 1) - 1, nrow = V)
  Y <- temp %*% errs
} else if(dist == "unif"){
  errs <-  t(chol(Omega.true)) %*% matrix(sample((-2):(2), size = n *V, replace = T)/ sqrt(2), nrow = V)
  Y <- temp %*% errs
}

B <- matrix(0, nrow = V, ncol = V)
Omega <- matrix(1, nrow = V, ncol = V)
tol <- 50

test.stat <- matrix(0, nrow = V, ncol = V)
for(i in 1:(V-1)){
  for(j in (i+1):V){
    
    B1 <- B
    B1[i, j] <- 1
    Omega[i, j] <- Omega[j, i] <- 0
    el_fit1 <- fitEL(Y = Y, B = B1, Omega = Omega,
                     meanEst = 1, type = "profile", high_moments = 4, optim.method = "BFGS")
    
    B2 <- B
    B2[j,i] <- 1
    el_fit2 <- fitEL(Y = Y, B = B2, Omega = Omega,
                     meanEst = 1, type = "profile", high_moments = 4, optim.method = "BFGS")
    
    diff <- el_fit2$lrt - el_fit1$lrt
    test.stat[i,j] <- diff
    cat("\n Testing: ")
    cat(paste(i, ",", j, "\n"))
#     if(i ==1 & j ==4){
#       mod1 <- el_fit1
#       mod2 <- el_fit2
#     }
    if( diff > tol){
      B <- B1
      cat(paste("ID Edge:", i, j ))
    } else if ( -diff > tol){
      B <- B2
      cat(paste("ID Edge:", j, i ))
    } else {
      cat("No Edge")
    }

  }
}

edit.dist[ind] <- sum((B - B.true.struct)^2)
}

B
Omega

B1 <- B
el_fit1 <- fitEL(Y = Y, B = B1, Omega = Omega,
                 meanEst = 1, type = "profile", high_moments = 4, optim.method = "BFGS")

B2 <- B
B2[4,1] <- 0
el_fit2 <- fitEL(Y = Y, B = B2, Omega = Omega,
                 meanEst = 1, type = "profile", high_moments = 4, optim.method = "BFGS")

el_fit3 <- fitEL(Y = Y, B = B.true.struct, Omega = Omega,
                 meanEst = 1, type = "profile", high_moments = 4, optim.method = "BFGS")

el_fit1$lrt
el_fit2$lrt