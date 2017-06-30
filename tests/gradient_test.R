#### Testing gradient for Naive EL ####
library(BCD)
library(numDeriv)



V <- 8
n <- 200
direct.num <- 8
bidirect.num <- 16
dist <- "t"

index <- c(1:V^2)
index <- index[lower.tri(matrix(0, nrow = V, ncol = V))]

B <- Omega <- matrix(0, nrow = V, ncol = V, byrow = T)
edge.list <- sample(index, size = direct.num + bidirect.num)
    
for(k in 1:direct.num){
  B[edge.list[k]] <- 1  
  }
if(bidirect.num > 0){
  for(k in c(1:bidirect.num)){
  Omega[edge.list[k +  direct.num]] <- 1
  }
}
Omega <- Omega + t(Omega) + diag(rep(1,V))
    
if(dist == "gauss"){
  B.true <- (rbinom(size = 1, n = V^2, prob = .5) *2 - 1) * runif(V^2, .2, 1) * B
  Omega.true <- matrix(0, nrow = V, ncol = V)
  Omega.true[(Omega == 1) & lower.tri(Omega)] <- rnorm(bidirect.num, 0, 1)
  Omega.true <- Omega.true + t(Omega.true)
  Omega.true <- Omega.true + diag(rep(1,V))
  diag(Omega.true) <- rowSums(abs(Omega.true)) + rchisq(V,1)
  
  temp <- solve(diag(rep(1,V)) - B.true)
  Y <- temp %*% t(MASS::mvrnorm(n = n, mu = rep(0, V), Sigma = Omega.true))
} else if(dist == "t") {
  t.df <- 5
  B.true <- (rbinom(size = 1, n = V^2, prob = .5) *2 - 1) * runif(V^2, .2, 1) * B
  Omega.true <- matrix(0, nrow = V, ncol = V)
  Omega.true[(Omega == 1) & lower.tri(Omega)] <- rnorm(bidirect.num, 0, 1)
  Omega.true <- Omega.true + t(Omega.true)
  Omega.true <- Omega.true + diag(rep(1,V))
  diag(Omega.true) <- rowSums(abs(Omega.true)) + rchisq(V,1)
  
  temp <- solve(diag(rep(1,V)) - B.true)
  Y <- temp %*% t(mvtnorm::rmvt(n = n, sigma = Omega.true * (t.df - 2) / t.df, df = t.df))
} else if(dist == "ln") {
  
  B.true <- (rbinom(size = 1, n = V^2, prob = .5) *2 - 1) * runif(V^2, .2, 1) * B
  sig.ln <- matrix(runif(V^2, -.1, .2), nrow = V) * Omega
  sig.ln <- sig.ln + t(sig.ln)
  diag(sig.ln) <- rowSums(abs(sig.ln)) + runif(V, .2, .5)
  errs <- t(exp(mvtnorm::rmvnorm(n, mean = rep(0, V), sigma = sig.ln)))
  errs <- errs - exp(.5 * diag(sig.ln))
  Omega.true <- exp(.5 * matrix((rep(diag(sig.ln), each = V) + rep(diag(sig.ln), V)), nrow = V)) * (exp(sig.ln) - 1)
  temp <- solve(diag(rep(1,V)) - B.true)
  
  Y <- temp %*% errs
}


out_ricf <- ricf(B, Omega, Y, maxIter = 1)
# 
vals <- c(B.true[which(B==1)], Omega.true[(Omega == 1) & lower.tri(Omega, diag = T)])
# vals <- c(out_naive_init$B[which(B==1)], out_naive_init$Omega[(Omega == 1) & lower.tri(Omega, diag = T)])
l_at_val <- sempl_input_naive(vals, Y, B, Omega, F, 1e-5, 100)
gr.analyt <- attr(l_at_val, "gradient")
gr.num <- grad(sempl_input_naive, vals, method = "Richardson", y_r =Y, b_r =  B, omega_r =  Omega,
     mean_est_r = F, tol =  1e-5, max_iter = 100)

sum((gr.analyt - gr.num)^2)

          