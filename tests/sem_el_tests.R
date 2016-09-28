#### Testing Empirical Likelihood for Structural Equation Models ####
library(microbenchmark)
set.seed(1000)
V <- 5
n <- 500

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

errs <- mvtnorm::rmvnorm(n, sigma = Omega.true)
Y <- solve(diag(rep(1, V)) - B.true) %*% t(errs)

microbenchmark(out_naive <- sempl(Y, B, Omega, meanEst = F, type = "naive"))
microbenchmark(out_prof <- sempl(Y, B, Omega, meanEst = F, type = "profile"))
