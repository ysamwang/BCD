#### BCD Test File ####

# True Values acyclic B
library(MASS)
library(microbenchmark)
set.seed(100)
V <- 5
n <- 10
B <- matrix(c(0,1,0,0,0,
              0,0,1,0,0,
              0,0,0,1,0,
              0,0,0,0,1,
              1,0,0,0,0), nrow = V, byrow = T)

Omega <- matrix(1, nrow = V, ncol = V)
Omega[1,3] <- Omega[3,1] <- 0
Omega[1,4] <- Omega[4,1] <- 0


B.true <- matrix(rnorm(V^2), nrow = V) * B
Omega.true <- ifelse(Omega,runif(V^2, -1, 1),0) + diag(rep(4, V))
Omega.true[lower.tri(Omega.true,diag=T)] = 0
Omega.true = Omega.true + t(Omega.true)
Omega.true = Omega.true + diag(rowSums(abs(Omega.true)) + rchisq(dim(Omega.true)[1], df = 1))

sigma = solve(diag(rep(1, V)) - B.true) %*% Omega.true %*% t(solve(diag(rep(1, V)) - B.true))
sum(eigen(sigma)$values < 0)

Y <- t(mvrnorm(n = n, mu = rep(0, V), Sigma = sigma))
# obj_given_B <- function(dual_var, B, Y, omega_rows, omega_cols)
# {
#   resid = (diag(rep(1, dim(B)[1])) - B) %*% Y
#   cross = resid[omega_rows, ] * resid[omega_cols, ]
#   pn = 1 / (dim(Y)[2]*(1 + dual_var %*% rbind(resid, cross) ))
#   return(-prod(pn))
# }
# 
# get_pn <- function(dual_var, B, Y, omega_rows, omega_cols)
# {
#   resid = (diag(rep(1, dim(B)[1])) - B) %*% Y
#   cross = resid[omega_rows, ] * resid[omega_cols, ]
#   pn = 1 / (dim(Y)[2]*(1 + dual_var %*% rbind(resid, cross) ))
#   return(pn)
# }
# 
# 
# 
# omega_rows <- c()
# omega_cols <- c()
# for(i in 2:V)
# {
#   for(j in 1:i)
#   {
#     if(Omega[i,j] == 0)
#     {
#       omega_rows <- c(omega_rows, i)
#       omega_cols <- c(omega_cols, j)
#     }
#   }
# }
# num_vars = V + sum(Omega[upper.tri(Omega, diag = F)] == 0)
# res <- optim(par = rep(0, num_vars), fn = obj_given_B, gr = NULL, B.true, Y, omega_rows, omega_cols, method = "BFGS")
# p_star <- get_pn(res$par, B.true, Y, omega_rows, omega_cols)
# 
# 
# 
# 

res <- sem_el_fitC(y_r = Y, b_r = B, omega_r = Omega,
            b_weights_r = B.true, d_r = rep(n, n),
            lambda_r = rep(0, V), gamma_r = rep(0, sum(Omega == 0)/2), v = V, tol = 1e-6, max_iter = 300)




