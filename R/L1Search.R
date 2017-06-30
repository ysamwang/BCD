sempl_L1 <- function(Y, lambdaVec, outerTol = 1e-4,
   innerTol = 1e-3, meanEst = F, maxIter = 100) {

  V <- dim(Y)[1]
  B.full <- matrix(1, nrow = V, ncol = V)
  diag(B.full) <- 0
  Omega <- diag(rep(1,V))

  ret <- list()

  out_ricf_init <- ricf(B = B.full, Omega = Omega, Y = Y - rowMeans(Y), BInit = NULL, OmegaInit = NULL, sigConv = 0, maxIter = 0, msgs = FALSE, omegaInitScale = .9)
  B.hat <- out_ricf_init$BHat

  for(l in 1:length(lambdaVec)){
    lambda <- lambdaVec[l]
    init_val <- c(B.hat[which(B.full == 1)])

    covarRestrict <- which(Omega == 0 & lower.tri(Omega), arr.ind = T) -1
    covarRestrict <- rbind(cbind(covarRestrict, rep(-1, dim(covarRestrict)[1]), rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,1], rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,2], rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,1], covarRestrict[,1]),
                           cbind(covarRestrict, covarRestrict[,2], covarRestrict[,2]))

    out_p <- lbfgs(call_eval = sempl_input_reg, call_grad = sempl_input_grad,
                 vars = init_val,y_r = Y, b_r = B.full,
                 mean_est_r = meanEst, covar_restrict = covarRestrict,
                 tol = innerTol, max_iter = maxIter, orthantwise_c = lambda,
                 epsilon = outerTol, invisible = 1)


    B.hat[which(B.full == 1)] <- out_p$par
    ret[[l]] <- list(lambda, B.hat)
  }
  return(ret)
}
