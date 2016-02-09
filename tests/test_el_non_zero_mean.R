  #### Testing Empirical Likelihood for Structural Equation Models ####
  
  set.seed(1000)
  V <- v<-  5
  
  B <- matrix(c(0, 0, 1, 0, 0,
                1, 0, 0, 0, 0,
                0, 1, 0, 0, 1,
                0, 0, 0, 0, 1,
                0, 0, 0, 0, 0) , nrow = V, ncol = V, byrow = T)
  Omega <- diag(rep(1,V))
  Omega[3,4] <- Omega[4,3] <- 1
  Omega[2,5] <- Omega[5,2] <- 1
  Omega.true <- diag(rep(1,V))
  Omega.true[3,4] <- Omega.true[4,3] <- .6
  Omega.true[2,5] <- Omega.true[5,2] <- -.2
  mu <- c(1, 2, 2, 1, 2)
  
  Y.true <- solve(diag(rep(1,V)) - B.true) %*% mu
  B.true <- runif(V^2, -1, 1) * B
  b.vec <- B.true[B == 1]
  o.vec <- Omega.true[lower.tri(Omega, diag = T) & Omega == 1]
  
  temp <- solve(diag(rep(1,V)) - B.true)
  sigma <- temp %*% Omega.true %*% t(temp)
  sigma
  
  sim.size <- 200
  n.list <- c(500)
  dist.list <- c("t")
  
  b.vec.ricf <- b.vec.el <- b.vec.el_no_mean <- matrix(0, nrow = sim.size, ncol = sum(B) * 3)
  o.vec.ricf <- o.vec.el <- o.vec.el_no_mean <- matrix(0, nrow = sim.size, ncol = (sum(Omega) - V)/2 + V)   
  mu.hat <- matrix(0, nrow = sim.size, ncol = V * 3)

  for(dist in dist.list){
    for(n.index in 1:length(n.list)){
      n <- n.list[n.index]
      cat(paste("\n",dist,n,":\n", sep = " "))
      for(sim in 1:sim.size)
      {
        cat(paste(sim, "."))
        
        #### Generate Data ####
        if(dist == "t"){
          t.df <- 5
          Y <- temp %*% (t(mvtnorm::rmvt(n, sigma = Omega.true * (t.df-2)/t.df, df = t.df)) + mu)
        } else if(dist == "lognormal"){
          sig.ln <- log(1 + (exp(1) - 1) * Omega.true)
          z <- t(MASS::mvrnorm(n = n, mu = rep(0, V), Sigma = sig.ln))
          Y <- temp %*% (exp(z) - exp(1/2))/ sqrt(exp(1) * (exp(1) - 1))
          
          } else if(dist == "gauss") {
          Y <- temp %*% t(MASS::mvrnorm(n = n, mu = mu, Sigma = Omega.true))
        }
        
        out_ricf <- ricf(B = B, Omega = Omega, Y = (Y - rowMeans(Y)), BInit = NULL,
                                     OmegaInit = NULL, sigConv = 0, maxIter = 1000,
                                     msgs = FALSE, omegaInitScale = .9)
        
        
        el_fit_null <- fitEL(Y = Y, B = B, Omega = Omega, meanEst = 1)
        el_fit_null_no_mean <- fitEL(Y = (Y - rowMeans(Y)), B = B, Omega = Omega, meanEst = 0)
  
        
        b.vec.ricf[sim,] <- out_ricf$BHat[B == 1]
        b.vec.el[sim,] <- el_fit_null$B.hat[B == 1]
        b.vec.el_no_mean[sim,] <- el_fit_null_no_mean$B.hat[B == 1]
  
  
        o.vec.ricf[sim,] <- out_ricf$OmegaHat[lower.tri(Omega, diag = T) & Omega == 1]
        
        o.vec.el[sim,] <- el_fit_null$Omega.hat[lower.tri(Omega, diag = T) & Omega == 1]
        o.vec.el_no_mean[sim,] <- el_fit_null_no_mean$Omega.hat[lower.tri(Omega, diag = T) & Omega == 1]
        mu.hat[sim, ] <- c(el_fit_null$mu.hat,
                           (diag(rep(1,V)) - el_fit_null_no_mean$B.hat) %*% rowMeans(Y),
                           (diag(rep(1,V)) - out_ricf$BHat) %*% rowMeans(Y))
                                                             
      }
    }
  #   saveRDS(error.b / sim.size, paste("record_b_",dist,".RDS", sep = ""))
  #   saveRDS(error.o / sim.size, paste("record_o_",dist,".RDS", sep = ""))
  #   saveRDS(reconstruct / sim.size, paste("record_sigma_",dist,".RDS", sep = ""))
  }
  
  
colMeans(recentered.mu <- t(t(mu.hat) - mu))
matrix(apply(recentered.mu, MAR = 2, FUN = sd), nrow = 3, byrow = T)

colMeans(t(t(o.vec.el) - o.vec)^2)
colMeans(t(t(o.vec.el_no_mean) - o.vec)^2)
colMeans(t(t(o.vec.ricf) - o.vec)^2)

matrix(apply(recentered.mu, MAR = 2, FUN = sd), nrow = 3, byrow = T)
