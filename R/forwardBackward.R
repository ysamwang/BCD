#' Discovery of Linear Structural Equation Models
#'
#' Produces estimated SEM from edge-wise tests
#'
#' @param Y data in V x N matrix
#' @param moments number of moments to test in estimating equations
#' @param fwdTol tolerance for difference on forward search
#' @param bwdTol tolerance for difference on backward search
#' @param fwdType type of fit for forward search: profile or euclid
#' @param bwdType type of fit for backward search: profile or euclid
#' @param meanEst whether or not to fit a mean
#' @param optimTol tolerance for optim
#' @param innerTol tolerance for dual problem
#' @param printStatus whether or not to print output
#' @export
forwardBackwardDiscovery <- function(Y, moments = 4, meanEst = F,
                                     outerTol = 1e-1, innerTol = 1e-2,
                                     printStatus = F, tieredSearch = F, vsNull = F,
                                     fwdTol = 30, bwdTol = 15){
  V <- dim(Y)[1]
  B <- matrix(0, nrow = V, ncol = V)

  # hold calculated differences
  lr.stat <- matrix(0, nrow = V, ncol = V)

  Omega <- diag(rep(1,2))
  B1 <- matrix(c(0, 1,
                 0, 0), nrow = 2, ncol = 2, byrow = T)
  B2 <- t(B1)
  B3 <- matrix(0, nrow = 2, ncol = 2)


  covarRestrict <- which(Omega == 0 & lower.tri(Omega), arr.ind = T) - 1
  if(moments == 3){
    covarRestrict <- rbind(cbind(covarRestrict, rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,1] ),
                           cbind(covarRestrict, covarRestrict[,2] ))
  } else if (moments == 4) {
    covarRestrict <- rbind(cbind(covarRestrict, rep(-1, dim(covarRestrict)[1]), rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,1] ,rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,2] ,rep(-1, dim(covarRestrict)[1])),
                           cbind(covarRestrict, covarRestrict[,1], covarRestrict[,1]),
                           cbind(covarRestrict, covarRestrict[, 2],covarRestrict[, 2] ))
  }

  #### Forward Search ####

    for(i in 2:V){
      for(j in 1:(i-1)){

        if(printStatus){
          cat(paste("\nForward:", i, j))
        }

        Y.test <- Y[c(i,j), ]






        ### Tests ####

        # edge i <- j
        el_fit1 <- sempl_cd(Y = Y.test, B = B1, Omega = Omega, meanEst = meanEst,
                            covarRestrict = covarRestrict, tol = innerTol,
                            outerTol = outerTol)

        # edge i -> j
        el_fit2 <- sempl_cd(Y = Y.test, B = B2,Omega = Omega, meanEst = meanEst,
                           covarRestrict = covarRestrict, tol = innerTol,
                           outerTol = outerTol)


        if(vsNull){
          # no edge i -\- j
          el_fit3 <- sempl_cd(Y = Y.test, B = B2,Omega = Omega, meanEst = meanEst,
                              covarRestrict = covarRestrict, tol = innerTol,
                              outerTol = outerTol)

          if (printStatus){
            cat(paste("\n\t",i, "," ,j, ": ", round(el_fit1,3), sep = ""))
            cat(paste("\n\t",j, "," ,i, ": ", round(el_fit2,3), sep = ""))
            cat(paste("\n\t",i, "," ,j, " (none) : ", round(el_fit3,3), sep = ""))
          }


          # update B if adding edge decreases LR stat by some tolerance
          if( (el_fit3 - min(el_fit2, el_fit1)) > fwdTol ){

            # add either edge i <- j or i -> j
            if(el_fit2 > el_fit1){
              B[i, j] <- 1
              lr.stat[i, j] <- (el_fit3 - min(el_fit2, el_fit1))
            } else {
              B[j, i] <- 1
              lr.stat[j, i] <- (el_fit3- min(el_fit2, el_fit1))
            }
          }


        } else {

          if (printStatus){
            cat(paste("\n\t",i, "," ,j, ": ", round(el_fit1,3), sep = ""))
            cat(paste("\n\t",j, "," ,i, ": ", round(el_fit2,3), sep = ""))
          }


          if (el_fit2 -  el_fit1  > fwdTol){
            B[i, j] <- 1
            lr.stat[i, j] <- el_fit2 -  el_fit1
          } else if (el_fit1 -  el_fit2  > fwdTol) {
            B[j, i] <- 1
            lr.stat[j, i] <- el_fit1 -  el_fit2
          }
      }

      }
    }



  order_of_prune <- cbind(which(lr.stat != 0 , arr.ind = T), lr.stat[which(lr.stat != 0)])
  order_of_prune <- order_of_prune[order(order_of_prune[,3]),]

  if(printStatus){
    cat("\nCurrent B:\n")
    print(order_of_prune)
    cat("\n")
  }


  Omega <- diag(rep(1, V))
  # Backward Pruning
  covarRestrict <- which(Omega == 0 & lower.tri(Omega), arr.ind = T) - 1
#   if(moments == 3){
#     covarRestrict <- rbind(cbind(covarRestrict, rep(-1, dim(covarRestrict)[1])),
#                            cbind(covarRestrict, covarRestrict[,1] ),
#                            cbind(covarRestrict, covarRestrict[,2] ))
#   } else if (moments == 4) {
#     covarRestrict <- rbind(cbind(covarRestrict, rep(-1, dim(covarRestrict)[1]), rep(-1, dim(covarRestrict)[1])),
#                            cbind(covarRestrict, covarRestrict[,1] ,rep(-1, dim(covarRestrict)[1])),
#                            cbind(covarRestrict, covarRestrict[,2] ,rep(-1, dim(covarRestrict)[1])),
#                            cbind(covarRestrict, covarRestrict[,c(1,1)]),
#                            cbind(covarRestrict, covarRestrict[,c(2,2)] ))
#   }

  edit <- T
  for(z in 1:length(order_of_prune[,1])){
        i <- order_of_prune[z, 1]
        j <- order_of_prune[z, 2]

        B1 <- B2 <- B
        B1[i, j] <- 0
        
        el_fit1 <- sempl_cd(Y = Y, B = B1, Omega = Omega, meanEst = meanEst,
                            covarRestrict = covarRestrict, tol = innerTol,
                            outerTol = outerTol)

      
        if(edit){
          el_fit2 <- sempl_cd(Y = Y, B = B2, Omega = Omega, meanEst = meanEst,
                              covarRestrict = covarRestrict, tol = innerTol,
                              outerTol = outerTol)
        }


        diff <- el_fit1 - el_fit2



            if(printStatus){
              cat(paste("Backward:", i, j))
              cat(paste(" ", round(diff,3), "\n"))
            }

            if( diff < bwdTol){
              B <- B1
              edit <- T
            } else {
              edit <- F
            }

          }

          return(B)
}
