# RICF with cycles v5.1
# Updated: April 9, 2014
# Fixed weird issue with norm that caused errors on old R installs

# INPUTS ---
# L, O: the 0-1 valued adjacency matrices indicating the structure of the graph G
# X: the n-by-p data matrix of observations from the true model
# Linit; Oinit: the initial values of the edge parameters for running the algorithm
# sigconv: Boolean (True = look for conv. in Sigma, False = look for conv. in L/O )
# tol: gives the max average entrywise absolute deviation in Sigma permited for convergence
# maxiter: gives the max number of iterations before divergence is accepted
# out: String: options: None/False, Final, All/True
# maxkap: maximum condition number accepted before error thrown
# B: (optional instead of L -- here B = t(L))
#
# OUTPUTS ---
# Sigmahat: the fitted value for Sigma resulting from the algorithm
# Bhat: the fitted value for the directed edge matrix (note: equal to I-t(Lambdahat))
# Omegahat, Lambdahat: the MLE for the parameter values resulting from the algorithm
# iterations: the number of iterations run by the algorithm before convergence or divergence accepted
# converged: TRUE or FALSE - based on whether or not the algorithm converged before maxiter reached
ricfR <- function(L = NULL, O, X, Linit = NULL, Oinit = NULL, sigconv=TRUE, tol=10^(-6),
                 maxiter=5000, out="none", maxkap = 10^9, B = NULL){
	if (!is.matrix(O) || !is.matrix(X))
		stop("O and X must be matrices!")
	if (nrow(O) != ncol(O))
		stop("O must be a square matrix!")
	if (sum(abs(O-t(O)))>0)
		stop("O must be a symmetric matrix!")
	if (!xor(is.null(L), is.null(B)))
		stop("Exactly one of B or L must be defined!")
	if (is.null(L)) {
		if (!is.matrix(B))
			stop("B must be a matrix!")
		if (nrow(B) != ncol(B))
			stop("B must be a square matrix!")
		L = t(B)
	}
	if (!is.matrix(L))
		stop("L must be a matrix!")
	if (nrow(L) != ncol(L) || nrow(O) != nrow(L))
		stop("L must be a square matrix of the same size as O!")
	p <- nrow(L)
	n <- nrow(X)
	# Initialize the directed edge parameters via OLS
	initL <- function(L, X) {
		Linit <- L
		parents <- apply(L, 2, function(x) which(x > 0))
		parcovariates <- lapply(parents, function(x) X[, x])
		len <- apply(L, 2, sum)
		for (j in 1:p) {
			if (len[j] > 0) {
				for (k in 1:len[j]) {
					mdl <- lm(as.numeric(X[,j])~as.numeric(X[, parents[[j]][k]]))
					Linit[parents[[j]][k], j] <- as.numeric(mdl$coef[2])
				}
			}
		}
		return(Linit)
	}
	# Initialize the bidirected edge parameters at random
	initO <- function(O) {
		R <- diag(p)
		for(i in 1:(p-1)) { for(j in (i+1):p) { R[j, i] <- R[i, j] <- rnorm(1) }}
		O2 <- R*O
		diag(O2) <- rep(0, p)
		diag(O2) <- rowSums(abs(O2)) + abs(rnorm(p))
		return(O2)
	}
	if (is.null(Linit)) { Linit <- initL(L, X) }
	if (is.null(Oinit)) { Oinit <- initO(O) }
	if (!is.matrix(Linit) || !is.matrix(Oinit))
		stop("O, X, Linit, and Oinit need to be matrices!")
	if (nrow(Linit) != ncol(Linit) || nrow(Oinit) != ncol(Oinit))
		stop("Linit and Oinit must be square matrices!")
	if (length(unique(c(nrow(Linit), nrow(Oinit), nrow(L), ncol(X)))) > 1)
		stop("One of the input matrices has the wrong dimension!")
	if (sum(abs(L * O)) != 0) {
		warning("The graph either contains a bow or a self-loop.",call. = FALSE, immediate. = TRUE)
	}	
	if (maxiter <= 0 || maxiter %% 1 != 0)
		stop("A positive integer is needed for the max number of iterations!")
	if (tol <= 0)
		stop("A positive tolerance is needed for convergence to be possible!")
	if (!is.logical(sigconv))
		stop("sigconv needs to take on a logical value!")
	out = tolower(out)
	if (!is.character(out))
		stop("Output needs to be a string: none/false, final, or all/true!")
	if (!(out == "true" || out == "all" || out == "final" || out == "false" || out == "none")) {
		stop("Output variable needs to be: none/false, final, or all/true!")
	}
	Det <- function(Lcur) {
		return(det(diag(nrow(Lcur)) - Lcur))
	}
	
	norm_vec <- function (x) sqrt(sum(x^2))
	Qp <- function(a) {
		p <- length(a)
		a2 <- a / norm_vec(a)
		Q <- qr.X(qr(a2), complete = TRUE)
		if(p == 1)
			return(Q)
		else
			return(cbind(Q[, 2:p], Q[, 1]))
	}
	Lcur <- Linit; Ocur <- Oinit
	iter <- 1
	repeat {
		for (i in 1:p) {
			pa <- which(L[, i] != 0)
			n.pa <- length(pa)
			sp <- which(O[, i] != 0)
			sp <- sp[sp != i]
			n.sp <- length(sp)
			len <- n.pa + n.sp
			IB <- diag(p) - Lcur
			Elessi <- (X %*% IB)[, -i]
			if (kappa(Ocur[-i, -i]) > maxkap) {
				stop(paste("The condition number of Ocur[-i, -i] is too large for node", i))
			}
			Zlessi <- Elessi %*% solve(Ocur[-i, -i])
			# The following line gets the indices of Zlessi corresponding to spouses
			zsp <- c(sp[sp < i], sp[sp > i] - 1)
			Zsp <- Zlessi[, zsp] # grab these columns in Zlessi
			Yi <- X[, i]
			Xi <- X[, pa]	
			if (len > 0) {
				if (n.pa > 0) {
					## DETERMINE a AND a0			
					a <- rep(0, n.pa)
					for (k in 1:n.pa) {
						temp <- Lcur
						temp[pa[k], i] <- 2
						det2 <- Det(temp)
						temp[pa[k], i] <- 1
						det1 <- Det(temp)
						a[k] <- det2 - det1
					}
					ind.pos <- which(a != 0)
					pa.pos <- pa[ind.pos]
					if (length(ind.pos) > 0) {
						
						## Build matrix Q and run the algorithm
						ind0 <- which(a == 0)
						pa0 <- pa[ind0]
						n.pa0 <- length(pa0)
						a <- a[ind.pos]
						temp[pa, i] <- 0
						a0 <- Det(temp)
						a0 <- a0 / sqrt(sum(a^2))
						a <- a / sqrt(sum(a^2))		
						Q <- Qp(a)
						M <- cbind(Zsp, X[, pa0], X[, pa.pos] %*% Q)
						qrM <- qr(M)
						Qtil <- t(qr.Q(qrM, complete = TRUE))
						R <- qr.R(qrM, complete = FALSE)	
						## Now solve for every ROTATED coef. (but NOT beta_n.pa)
						coef <- rep(0, len)
						if (len > 1) {coef[-len] <- (Qtil %*% Yi)[1:(len - 1)]}
						## Now solve for the ROTATED beta_n.pa
						y0 <- sum(((Qtil %*% Yi)^2)[(len + 1):n])
						r <- R[len, len]
						yp <- (Qtil %*% Yi)[len]
            sol <- (y0 + a0 * yp * r + yp^2) / (r * (a0 * r + yp))
           
						coef[len] <- r * sol
						
						## Convert the coefficients back (multiply by R inverse)
						coef.new <- solve(R, coef)
						## Multiply by Q to find betas optimizing the origninal function
						Ocur[i, sp] <- Ocur[sp, i] <- if(n.sp > 0){coef.new[1:n.sp]} else{c()}
						Lcur[pa0, i] <- if(n.pa0 > 0){coef.new[(n.sp + 1):(n.sp + n.pa0)]} else{c()}
						Lcur[pa.pos, i] <- Q %*% (coef.new[(n.sp + n.pa0 + 1):len])		
					}
					else {
						## Perform linear regression for Betas and (possibly) Omegas
						mdl <- lm(Yi ~ cbind(Zsp, Xi) - 1)
						coef <- mdl$coef
						if (n.sp > 0) {
							if (any(is.na(coef))) {
								stop(paste("Collinearity observed for node", i))
							}
							Ocur[sp, i] <- Ocur[i, sp] <- coef[1:n.sp]
						}
						Lcur[pa, i] <- coef[(n.sp + 1):len]
					}
				}
				else {
					## Solve a simplified version with no Betas (Only Omegas)
					mdl <- lm(Yi ~ Zsp - 1)
					Ocur[sp, i] <- Ocur[i, sp] <- mdl$coef
				}
			}	
			## Find the variance omega_{ii}
			res <- (Yi - as.matrix(Xi) %*% Lcur[pa, i] - as.matrix(Zsp) %*% Ocur[sp, i])
			RSS <- t(res) %*% res
			if (kappa(Ocur[-i, -i]) > maxkap) {
				stop(paste("The Condition number of Ocur[-i, -i] is too large for node", i))
			}
			
			### print ###
			
			
			
			Ocur[i,i] <- (RSS/n) + Ocur[i, -i] %*% solve(Ocur[-i, -i]) %*% Ocur[-i, i]
			
		}	
		sigcur <- t(solve(diag(p) - Lcur)) %*% Ocur %*% solve(diag(p) - Lcur)
		bhat <- diag(p)-t(Lcur)
		if (iter == 1) {
			if (out == "true" || out == "all"){cat(iter, "\n")}
			if (maxiter == 1) {
				if (out == "true" || out == "all" || out == "final") {
					cat("Sigmahat \n")
					print(sigcur)
					cat("\nBhat \n")
					print(bhat)
					cat("\nOmegahat \n")
					print(Ocur)
					cat("\nLambdahat \n")
					print(Lcur)
					cat("\niterations \n")
					print(iter)
				}
				break
			}
		}
		else if (iter > 1) {
			dsig <- mean(abs(sigcur - sigpast))
			dLO <- sum(abs(Lcur - Lpast) + abs(Ocur - Opast)) / (sum(L) + sum(O))
			if (out == "true" || out == "all") {
				dsig6 <- format(round(dsig, 6), nsmall = 6)
				dLO6 <- format(round(dLO, 6), nsmall = 6)
				cat(iter, " Avg Change in L & O: ", dLO6, "| Avg Change in Sigma: ", dsig6, "\n")
			}

			if ((sigconv && dsig < tol) || (!sigconv && dLO < tol) || iter >= maxiter) {
				if (out == "true" || out == "all" || out == "final") {
					cat("Sigmahat \n")
					print(sigcur)
					cat("\nBhat \n")
					print(bhat)
					cat("\nOmegahat \n")
					print(Ocur)
					cat("\nLambdahat \n")
					print(Lcur)
					cat("\niterations \n")
					print(iter)					
				}
				break
			}
		}
		sigpast <- sigcur; Lpast <- Lcur; Opast <- Ocur
		iter <- iter + 1
	}
	return(list(Sigmahat = sigcur, Bhat = bhat, Omegahat = Ocur, Lambdahat = Lcur, iterations = iter, converged = (iter < maxiter)))
}




# Optional method for calling ricf, using a mag object instead of L and O.
# Note: undirected edges not allowed here.
ricfmag <- function(mag, X, Linit = NULL, Oinit = NULL, sigconv=TRUE, tol=10^(-6), maxiter=5000, out="none", maxkap = 10^9)
{
	if (!is.matrix(mag) || !is.matrix(X))
		stop("mag and X must be matrices!")
	if (nrow(mag) != ncol(mag))
		stop("mag must be a square matrix!")	
	L <- mag %% 100
	O <- mag %/% 100 + diag(nrow(mag))
	return(ricf(L,O,X,Linit,Oinit,sigconv,tol,maxiter,out,maxkap))
}