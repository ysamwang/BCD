# returns covar argument for sem::specifyModel function
# Takes Omega Matrix of 1/0 specifying structural 0's in Bidirected Graph
getCovar <- function(Omega){
  p <- nrow(Omega)
  covar <- c()
  for(i in 1:p)
  {
    for(j in 1:i){
      if(i == j){
        covar <- c(covar, paste("x", i, sep = ""))
      } else {
        if(Omega[i,j] ==1){
          covar <- c(covar, paste(paste("x", i, sep =""), paste("x", j, sep=""), sep = ", "))
        }
      }
      
    }
  }
  return(covar)
}

# returns text argument for sem::specifyModel function
# Takes B Matrix of 1/0 specifying structural 0's in Directed Graph
# B_ij = 1 implies a directed edge from j -> i
getMod <- function(B){
  p <- nrow(B)
  mod <- ""
  first <- 1
  for(i in 1:p){
    for(j in 1:p){
      if(B[i,j] == 1){
        if(first){
          mod <- paste("x", j, "->", "x",i,", beta",i,j, sep = "")
          first <- 0 
        } else {
          mod <- paste(mod, "\n ","x", j, "->", "x",i,", beta",i,j, sep = "")
        }
      }
    }
  }
  return(mod)
}

specifyModelSEM <- function(B, Omega){
  return(sem::specifyModel(text = getMod(B), exog.variances = T, endog.variances = T, covs = getCovar(Omega), quiet = T))
}

microbenchmark::microbenchmark(out.sem <- sem(model = mod, data = dat, objective = objectiveGLS))