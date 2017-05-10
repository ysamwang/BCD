
tau.stat <- function(B, A, k){
  mean(A^(k - 1) * B) * mean(A^2) - mean(A^k) * mean(A * B)
}

dag_search_non_gaussian <- function(Y, K = 3, cutoff){
  ordering <- c()
  V <- dim(Y)[2]
  
  Y <- t(t(Y) - rowMeans(Y))
  
  remaining.edges <- c(1:V)
  
  # Get Ordering
  while(length(remaining.edges) > 1) {
    tau.min <- rep(0, length(remaining.edges))
    
    for(i in 1:length(remaining.edges)){
        A <- Y[, remaining.edges[i]]
        B <- Y[, remaining.edges[-i]]
        tau.min[i] <- min(apply(B, MAR = 2, tau.stat, A, K)) 
    }
    root <- which.min(abs(tau.min))
    ordering <- c(ordering, remaining.edges[root])
    remaining.edges <- remaining.edges[-root]
  }
  return(ordering)
  
  
}




