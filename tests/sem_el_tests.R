#### Testing Empirical Likelihood for Structural Equation Models ####
library(microbenchmark)
set.seed(1000)
V <- 5

B <- matrix(c(0, 0, 0, 0, 0,
              0, 0, 0, 0, 0,
              1, 1, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 1, 0, 0)
            , nrow = V, ncol = V, byrow = T)

Omega <- diag(rep(1,V))
Omega[1, 2] <- Omega[2, 1] <- 1
Omega[1, 4] <- Omega[4, 1] <- 1
Omega[5, 2] <- Omega[2, 5] <- 1

Omega.true <- diag(rep(1, V))
Omega.true[1, 2] <- Omega.true[2, 1] <- .6
Omega.true[1, 4] <- Omega.true[4, 1] <- -.2
Omega.true[5, 2] <- Omega.true[2, 5] <- .4

B.true <- round(runif(V^2, -1, 1),1) * B



sem.mod <- specifyModelSEM(B, Omega)

b.vec <- B.true[B == 1]
o.vec <- Omega.true[lower.tri(Omega, diag = T) & Omega == 1]

temp <- solve(diag(rep(1,V)) - B.true)
sigma <- temp %*% Omega.true %*% t(temp)
sigma

sim.size <- 500
n.list <- c(20, 30, 50, 100, 500)
dist.list <- c("gauss", "t", "lognormal")


for(dist in dist.list){
  error.b <- matrix(0, ncol = sum(B) * 3, nrow = length(n.list))
  error.o <- matrix(0, ncol =(V + (sum(Omega == 1) - V)/2) * 3, nrow = length(n.list))
  converged <- time <-  matrix(0, ncol = 3, nrow = length(n.list))
  for(n.index in 1:length(n.list)){
    n <- n.list[n.index]
    cat(paste("\n",dist,n,":\n", sep = " "))
    for(sim in 1:sim.size)
    {
      cat(paste(sim, "."))
      
      #### Generate Data ####
      if(dist == "t"){
        t.df <- 5
        Y <- temp %*% (t(mvtnorm::rmvt(n, sigma = Omega.true * (t.df-2)/t.df, df = t.df)))
      } else if(dist == "lognormal"){
        sig.ln <- log(1 + (exp(1) - 1) * Omega.true)
        z <- t(MASS::mvrnorm(n = n, mu = rep(0, V), Sigma = sig.ln))
        Y <- temp %*% ((exp(z) - exp(1/2))/ sqrt(exp(1) * (exp(1) - 1)))
      } else if(dist == "gauss") {
        Y <- temp %*% (t(MASS::mvrnorm(n = n, mu = rep(0, V), Sigma = Omega.true)))
      }
      
      time.el <- microbenchmark(out_el <- fitEL(Y = Y, B = B, Omega = Omega,
                                                       meanEst = 0, naive = FALSE),
                               times = 2, control = list(warmup = 2))
      
      
      time.ricf <- microbenchmark(out_ricf <- ricf(B = B, Omega = Omega, Y = Y, maxIter = 1000,
                                                 omegaInitScale = .9),
                                times = 2, control = list(warmup = 2))
      
      dat <- as.data.frame(t(Y))
      names(dat) <- paste("x", c(1:V), sep = "")

      time.sem <- microbenchmark(out_sem <- sem(model = sem.mod, data = dat),
                                 times = 2, control = list(warmup = 2))
       

      b.vec.el <- out_el$B.hat[B == 1]
      o.vec.el <- out_el$Omega.hat[lower.tri(Omega, diag = T) & Omega == 1]
      
      b.vec.ricf <- out_ricf$BHat[B == 1]
      o.vec.ricf <- out_ricf$OmegaHat[lower.tri(Omega, diag = T) & Omega == 1]
      
      b.vec.gls <- out_sem$A[B == 1]
      o.vec.gls <- out_sem$P[lower.tri(Omega, diag = T) & Omega == 1]
      
      
      error.b[n.index, ] <- error.b[n.index, ] + c((b.vec.el - b.vec)^2 * el_fit_reduced$converged,
                                                   (b.vec.ricf - b.vec)^2 * el_fit_naive$converged,
                                                   (b.vec.gls - b.vec)^2 * el_fit_naive$converged)
      
      error.o[n.index, ] <- error.o[n.index, ] + c((o.vec.el - o.vec)^2 * el_fit_reduced$converged,
                                                   (o.vec.ricf - o.vec)^2 * el_fit_naive$converged,
                                                   (o.vec.gls - o.vec)^2 * el_fit_naive$converged)
      
      converged[n.index, ] <- converged[n.index, ] + c(el_fit_reduced$converged,
                                                       el_fit_naive$converged)
      time[n.index, ] <- time[n.index, ] + c(mean(time.r$time) * el_fit_reduced$converged,
                                             mean(time.n$time) * el_fit_naive$converged)
      
    }
  }
  saveRDS(error.b, paste("record_b_",dist,".RDS", sep = ""))
  saveRDS(error.o, paste("record_o_",dist,".RDS", sep = ""))
  saveRDS(converged, paste("record_converged_",dist,".RDS", sep = ""))
  saveRDS(time, paste("record_time_",dist,".RDS", sep = ""))
}


