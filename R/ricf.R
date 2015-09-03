#' Fitting Cyclic Linear Structural Equation Models
#' 
#' 
#' Estimates MLE's for cyclic linear SEM's as described in DFW
#' 
#'  
#' @param B V by V matrix with {0,1} giving structure of directed edges
#' @param Omega V by V matrix with {0,1} giving structure of bi-directed edges
#' @param BInit V by V matrix giving initial edges weights for directed edges
#' @param OmegaInit V by V matrix giving initial edge weights for bi-directed edges
#' @param signConv boolean which specifies how to measure convergence. \code{TRUE} looks for
#'    convergence in Sigma while \code{FALSE} looks for convergence in the actual edge weight estimates  
#' @param tol convegence tolerance
#' @param maxIter integer specifying the maximum number of iterations
#' @return \item{sigmaHat}{estimated covariance matrix at convergence}
#'    \item{bHat}{estimated B matrix (edge weights for directed edges) at converegence}
#'    \item{omegaHat}{estimated Omega (edge weights for bi-directed edges) at convergence}
#'    \item{iterations}{number of iterations until convergence. a single iteration is considered
#'    a pass through all nodes}
#'    \item{converged}{boolean on whether or not the algorithm converged before the max iterations}
ricf <- function(B, Omega, Y, BInit = NULL , OmegaInit = NULL, sigConv = TRUE,
                 tol=10e-6, maxIter=5000) {
  
  if (!is.matrix(Omega) || !is.matrix(Y) || !is.matrix(B))
    stop("Omega B, and X must be matrices!")
  if (nrow(Omega) != ncol(Omega))
    stop("Omega must be a square matrix!")
  if (sum(abs(Omega-t(Omega)))> 0)
    stop("O must be a symmetric matrix!")
  if (nrow(B) != ncol(B) || nrow(Omega) != nrow(B))
      stop("B must be a square matrix! the same size as Omega")
  
  if (sum(abs(B * Omega)) != 0) {
    warning("The graph either contains a bow or a self-loop.",call. = FALSE, immediate. = TRUE)
  }	
  if (maxIter <= 0 || ceiling(maxIter) != maxIter)
    stop("A positive integer is needed for the max number of iterations!")
  if (tol <= 0)
    stop("A positive tolerance is needed for convergence to be possible!")

  ret <- bcdC(B, Omega, BInit, OmegaInit, Y, maxIter = maxIter, sigConv = F, maxKap = 1e9, tol= 1e-6)
  
  return(ret)

}
