set.seed(101)

V <- 3
B <- matrix(c(0,1,0,
              0,0,1,
              0,0,0), nrow = V, ncol = V, byrow = T)
Omega <- diag(rep(1, V))

B1 <- B
Omega1 <- Omega

B2 <- matrix(c(0, 1, 0,
               0, 0, 0,
               0, 0, 0), byrow = T, nrow = V)
Omega2 <- diag(rep(1, V))
Omega2[2,3] <- Omega2[3,2] <- 1



B3 <- matrix(c(0, 0, 0,
               0, 0, 0,
               0, 0, 0), byrow = T, nrow = V)
Omega3 <- Omega
Omega3[1,2] <- Omega3[2,1] <- 1
Omega3[2,3] <- Omega3[3,2] <- 1

B4 <- matrix(c(0, 0, 0,
               1, 0, 0,
               0, 0, 0), byrow = T, nrow = V)
Omega4 <- diag(rep(1, V))
Omega4[2,3] <- Omega4[3,2] <- 1

n <- 1000

sim.size <- 100

B.true <- round(rnorm(V^2),2) * B
Omega.true <- diag(rep(1, V))
# Omega.true[1,2] <- Omega.true[2,1] <- .5
# Omega.true[3,2] <- Omega.true[2,3] <- -.2

temp <- solve(diag(rep(1,V)) - B.true)
sig.ln <- log(1 + (exp(1) - 1) * Omega.true)

lrt1 <- lrt2 <- lrt3 <- lrt4 <-  rep(0, sim.size)
mse1.b <- mse2.b <- mse3.b <- rep(0, sim.size)
mse1.o <- mse2.o <- mse3.o <- rep(0, sim.size)

sigma <- temp %*% Omega.true %*% t(temp)
dist <- "lognormal"


library(distr)
D <- DExp(rate = sqrt(2))
for(i in 1:sim.size){
  cat("\n")
  cat(i)
  
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
    errs <-  t(chol(Omega.true)) %*% matrix(r(D)(n * V), nrow = V)
    Y <- temp %*% errs
  } else if(dist == "gamma"){
    errs <-  t(chol(Omega.true)) %*% matrix(rgamma(n * V, 1, 1) - 1, nrow = V)
    Y <- temp %*% errs
  } else if(dist == "unif"){
    errs <-  t(chol(Omega.true)) %*% matrix(sample((-2):(2), size = n *V, replace = T)/ sqrt(2), nrow = V)
    Y <- temp %*% errs
  }
  
  el_fit1 <- fitEL(Y = Y, B = B1, Omega = Omega1,
                   meanEst = 0, type = "profile", high_moments = 3, optim.method = "BFGS")
  
  
  
  el_fit2 <- fitEL(Y = Y, B = B2, Omega = Omega2,
                   meanEst = 0, type = "profile", high_moments = 3, optim.method = "BFGS")
  
  el_fit3 <- fitEL(Y = Y, B = B3, Omega = Omega3,
                   meanEst = 0, type = "profile", high_moments = 3, optim.method = "BFGS")

  el_fit4 <- fitEL(Y = Y, B = B4, Omega = Omega4,
                   meanEst = 0, type = "profile", high_moments = 3, optim.method = "BFGS")
  
#   el_fit3 <- fitEL(Y = Y, B = B3, Omega = Omega,
#                    meanEst = 0, naive = FALSE, high_moments = high_moments, maxInnerIter = 1000)
#   
#   el_fit4 <- fitEL(Y = Y, B = B4, Omega = Omega4,
#                    meanEst = 0, naive = FALSE, high_moments = high_moments, maxInnerIter = 1000)
#   
#   el_fit5 <- fitEL(Y = Y, B = B5, Omega = Omega5,
#                    meanEst = 0, naive = FALSE, high_moments = high_moments, maxInnerIter = 1000)
  
  lrt1[i] <- el_fit1$lrt * (abs(sum(el_fit1$p) - 1) < 1e-10)
  lrt2[i] <- el_fit2$lrt * (abs(sum(el_fit2$p) - 1) < 1e-10)
  lrt3[i] <- el_fit3$lrt * (abs(sum(el_fit3$p) - 1) < 1e-10)
  lrt4[i] <- el_fit4$lrt * (abs(sum(el_fit4$p) - 1) < 1e-10)
 
#   mse1.b[i] <- mean((el_fit1$B.hat - B.true)^2)
#   mse2.b[i] <- mean((el_fit2$B.hat - B.true)^2)
#   mse3.b[i] <- mean((el_fit3$B.hat - B.true)^2)
#   
#   mse1.o[i] <- mean((el_fit1$Omega.hat - Omega.true)^2)
#   mse2.o[i] <- mean((el_fit2$Omega.hat - Omega.true)^2)
#   mse3.o[i] <- mean((el_fit3$Omega.hat - Omega.true)^2)
  
}


# png("causal_discovery_9.png", width = 700, height = 500)
par(mfrow = c(1,4), oma = c(0,0,2,0))
hist(lrt1, main = "LR Stat")
mtext("True Model")

hist(lrt2 - lrt1, main = "LR Difference", xlab = "1<-2<-3")
mtext(paste(round(mean((lrt2 - lrt1) >= 0), 2), "Correct"))

hist(lrt3 - lrt1, main = "LR Difference", xlab = "1<-2->3")
mtext(paste(round(mean((lrt3 - lrt1) >= 0), 2), "Correct"))

hist(lrt3 - lrt1, main = "LR Difference", xlab = "1<-2->3")
mtext(paste(round(mean((lrt4 - lrt1) >= 0), 2), "Correct"))

mtext("True Model: 1 -> 2 -> 3; n = 1000; T distribution (4th moment)", outer = T)
# dev.off()

