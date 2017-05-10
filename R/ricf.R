#' Fitting Cyclic Linear Structural Equation Models
#' 
#' 
#' Estimates MLE's for cyclic linear SEM's as described in DFW
#' 
#'  
#' @param B V by V matrix with {0,1} giving structure of directed edges. B[i, j] = 1 indicates the edge j -> i.
#' @param Omega V by V matrix with {0,1} giving structure of bi-directed edges
#' @param BInit V by V matrix giving initial edges weights for directed edges. If BInit is NULL,
#'    a default initialization will be used. 
#' @param OmegaInit V by V matrix giving initial edge weights for bi-directed edges
#' @param signConv boolean which specifies how to measure convergence. \code{TRUE} looks for
#'    convergence in Sigma while \code{FALSE} looks for convergence in the actual edge weight estimates  
#' @param tol convegence tolerance
#' @param maxIter integer specifying the maximum number of iterations
#' @param msgs boolean on whether to print warning messages to command line if there are bows in the graph
#' @param omegaInitScale value in (0,1) which determines how to initialize estimates for the bidirected edges.
#'    The directed edges are initialized through OLS regression and the bidirected edge weights are initialized
#'    by placing the structural 0's in the covariance of the residuals. If the resulting initialization is not positive
#'    definite, we scale the elements of corresponding row/columns such that sum(OmegaInit[i,-i]) =  OmegaInit[i,i] * omegaInitScale
#' @return \item{sigmaHat}{estimated covariance matrix at convergence}
#'    \item{bHat}{estimated B matrix (edge weights for directed edges) at converegence}
#'    \item{omegaHat}{estimated Omega (edge weights for bi-directed edges) at convergence}
#'    \item{iterations}{number of iterations until convergence. a single iteration is considered
#'    a pass through all nodes}
#'    \item{converged}{boolean on whether or not the algorithm converged before the max iterations}
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
