#' Fitting Cyclic Linear Structural Equation Models
#' 
#' 
#' Estimates MLE's for cyclic linear SEM's as described in Drton, Fox, Wang (2019)
#' 
#'  
#' @param B V by V matrix with {0, 1} giving structure of directed edges. B[i, j] = 1 indicates the edge j -> i.
#' @param Omega V by V matrix with {0, 1} giving structure of bi-directed edges. Omega[i,h] = 1 indicates the edge i <-> j.
#'    The diagonal elements should also be 1
#' @param Y V by n data matrix where each row corresponds to an observed variable and each column 
#'    corresponds to a multivariate observation. The method assumes that each variable (row) is mean 0.  
#' @param BInit V by V matrix giving initial edges weights for directed edges. If BInit is NULL,
#'    a default initialization will be used. 
#' @param OmegaInit V by V matrix giving initial edge weights for bi-directed edges. If OmegaInit 
#'    is NULL, a default initialization will be used
#' @param sigConv boolean which specifies how to measure convergence. \code{TRUE} checks for
#'    convergence in Sigma while \code{FALSE} checks for convergence in the edge weight estimates  
#' @param tol convegence tolerance
#' @param maxIter integer specifying the maximum number of iterations
#' @param msgs boolean on whether to print warning messages to command line if there are bows in the graph
#' @param omegaInitScale value in (0,1) which determines how to initialize estimates for the bidirected edges.
#'    The directed edges are initialized through OLS regression and the bidirected edge weights are initialized
#'    by placing the structural 0's in the covariance of the residuals. If the resulting initialization is not positive
#'    definite, we scale the elements of corresponding row/columns such that sum(OmegaInit[i,-i]) =  OmegaInit[i,i] * omegaInitScale
#' @return \item{SigmaHat}{estimated covariance matrix at convergence}
#'    \item{OmegaHat}{estimated Omega (edge weights for bi-directed edges) at convergence}
#'    \item{BHat}{estimated B matrix (edge weights for directed edges) at converegence}
#'    \item{Iter}{number of iterations until convergence. a single iteration is considered
#'    a pass through all nodes}
#'    \item{Converged}{boolean on whether or not the algorithm converged before the max iterations}
#'    
#' @examples
#' ## Select True Parameters
# Omega.weights <- diag(rep(1, 4))
# Omega.weights[1, 2] <- Omega.weights[2, 1] <- .5
# B.weights <- matrix(0, nrow = 4, ncol = 4)
# B.weights[3, 1] <- .5; B.weights[4, 2] <- .7
# ## Generate data
# n <- 200
# epsilon <- t(mvtnorm::rmvnorm(n, mean = rep(0, 4), sigma = Omega.weights))
# Y <- solve(diag(rep(1, 4)) - B.weights, epsilon)
#
# ## Form ricf arguments
#' B <- (B.weights != 0) + 0
#' Omega <- (Omega.weights != 0) + 0 
#' ricf(B, Omega, Y)
#' out <- ricf(B, Omega, Y)
#' se <- var.ricf(Y, B, Omega,  out$BHat, out$OmegaHat, type = "expected")
#' # get parameter estimates and marginal standard errors
#' output <- cbind(c(out$BHat[which(B !=0)], out$OmegaHat[which(Omega !=0 & lower.tri(Omega, diag = T))]),
#'                sqrt(diag(se)))
#' colnames(output) <- c("estimates", "SE")
#' rownames(output) <- rownames(se)
#' output

ricf <- function(B, Omega, Y, BInit = NULL , OmegaInit = NULL, sigConv = TRUE,
                 tol=1e-7, maxIter=5000, msgs = TRUE, omegaInitScale = .9, maxKap = 1e13) {
  
  if (!is.matrix(Omega) || !is.matrix(Y) || !is.matrix(B))
    stop("Omega B, and X must be matrices!")
  if (nrow(Omega) != ncol(Omega))
    stop("Omega must be a square matrix!")
  if (sum(abs(Omega-t(Omega)))> 0)
    stop("O must be a symmetric matrix!")
  if (nrow(B) != ncol(B) || nrow(Omega) != nrow(B))
      stop("B must be a square matrix! the same size as Omega")
  
  if (msgs & sum(abs(B * Omega)) != 0) {
    # only print if msgs = TRUE
    warning("The graph either contains a bow or a self-loop.",call. = FALSE, immediate. = TRUE)
  }	
  if (maxIter < 0 || ceiling(maxIter) != maxIter)
    stop("A positive integer is needed for the max number of iterations!")
  if (tol <= 0)
    stop("A positive tolerance is needed for convergence to be possible!")

  # If All checks run, then call bcdC function
  ret <- bcdC(Br = B, Omegar = Omega, BInitr = BInit, OmegaInitr = OmegaInit, Yr = Y,
              maxIter = maxIter, sigConv = sigConv, maxKap = maxKap, tol = tol, omegaInitScale = omegaInitScale)
  
  return(ret)

}
