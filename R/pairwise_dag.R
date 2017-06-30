sempl_dag_search <- function(Y, moments = 3, cutoff = 3){
  V <- dim(Y)[1]
  B <- matrix(0, nrow = V, ncol = V)
  
  ordered <- c()
  unordered <- c(1:V)
  lrt <- rep(0, length(unordered))
  Y.inner <- Y
  
  while(length(unordered) > 1) {
  
    Y.inner <- Y - B %*% Y
    Y.inner <- Y.inner[unordered, ]
    for(k in 1:length(unordered)){
      print(k)
      B.inner <- matrix(0, nrow = length(unordered), ncol = length(unordered))
      B.inner[, k] <- 1
      B.inner[k, k] <- 0
      Omega.inner <- matrix(1, nrow = length(unordered), ncol = length(unordered))
      Omega.inner[, k] <- Omega.inner[k, ] <- 0
      Omega.inner[k, k] <- 1 
      
      covarRestrict <- .getCovarRestrict(Omega.inner, moments = moments)
      out <- sempl(Y.inner, B.inner, Omega.inner, meanEst = F,
                   covarRestrict = covarRestrict)
      lrt[k] <- out$lr
      names(out)
    }
    
  
    k <- which.min(lrt)
    B.inner <- matrix(0, nrow = length(unordered), ncol = length(unordered))
    B.inner[, k] <- 1
    B.inner[k, k] <- 0
    Omega.inner <- matrix(1, nrow = length(unordered), ncol = length(unordered))
    Omega.inner[, k] <- Omega.inner[k, ] <- 0
    Omega.inner[k, k] <- 1 
    
    covarRestrict <- .getCovarRestrict(Omega.inner, moments = moments)
    out <- sempl(Y.inner, B.inner, Omega.inner, meanEst = F,
                 covarRestrict = covarRestrict)
    B[unordered, unordered[k]] <- out$B.hat[, k]
    
    ordered <- c(ordered, unordered[k])
    unordered <- unordered[-k]
    
  }
  
  B.trial <- ifelse(abs(B) > 0, 1, 0)
  Omega.trial <- diag(rep(1,V))
  out <- sempl(Y, B.trial, Omega.trial, meanEst = F)
  v <- var.el(Y, B.trial, Omega.trial, out$B, out$Om, p = out$p)
  
  t.val <- out$B[which(B.trial == 1)] / sqrt(diag(v)[1:sum(B.trial)])
  sig <- which(abs(t.val) > cutoff)
  
  B <- matrix(0, nrow = V, ncol = V)
  B[which(B.trial == 1, arr.ind = T)[sig,]] <- 1
  return(B)
}

getCovarRestrict <- function(Omega, moments = 3){
  covar <- which(Omega == 0 & lower.tri(Omega), arr.ind = T) - 1
  ret <- cbind(covar, matrix(-1, nrow = dim(covar)[1], ncol = moments - 2))
  for(i in 3:moments){
    ret <- rbind(ret,
                 cbind(covar, matrix(covar[, 1], nrow = dim(covar), ncol = i - 2),
                       matrix(-1, nrow = dim(covar), ncol = moments - i)),
                 cbind(covar, matrix(covar[, 2], nrow = dim(covar), ncol = i - 2),
                       matrix(-1, nrow = dim(covar), ncol = moments - i))
                 )
  }
  return(ret)
}  
  
  