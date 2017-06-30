
n.list <- c(100, 250, 500, 1000)
dist.list <- c("gauss", "t", "poisson", "lognormal")
type <- dist.list[4]

b.error <- omega.error <- matrix(0,nrow = 3, ncol = length(n.list))
 
for(n in 1:length(n.list))
{
  rec <- readRDS(paste("record_", type, "_", n.list[n],".RDS", sep = ""))
  
  valid <- (rec$times > 0) & (rec$times_euclid > 0 ) & (rec$times_ricf >0)
  
  b.error[,n] <-c(median(rec$error_b[valid]) ,
                  median(rec$error_b_ricf[valid]),
                  median(rec$error_b_euclid[valid]))
  
  omega.error[,n] <- c(median(rec$error_omega[valid]) ,
                       median(rec$error_omega_ricf[valid]),
                       median(rec$error_omega_euclid[valid]))
}

y.max <- ifelse(max(b.error[,1]) > min(b.error[,1])*5, max(b.error[,2]) * 1.5, max(b.error))
plot(-1, -1, xlim = c(0,5), ylim = c(0, y.max), xaxt = "n",
     ylab = "Mean Squared Error", main = paste(type, "Distributed Errors MSE"),
     xlab = "Sample Size")
mtext("Directed Edge Weights")
points(b.error[1,], pch = "E")
points(b.error[2,], pch = "R")
points(b.error[3,], pch = "U")
axis(1, labels = c(100,250, 500,1000), at = c(1:4))

y.max <- ifelse(max(omega.error[,1]) > min(omega.error[,1])*5, max(omega.error[,2])*1.5,
                max(omega.error))
plot(-1, -1, xlim = c(0,5), ylim = c(0, y.max), xaxt = "n",
     ylab = "Median Squared Error", main = paste(type, "Distributed Errors MSE"),
     xlab = "Sample Size")
mtext("Bi-directed Edge Weights")
points(omega.error[1,], pch = "E")
points(omega.error[2,], pch = "R")
points(omega.error[3,], pch = "U")
axis(1, labels = c(100,250, 500,1000), at = c(1:4))

library(xtable)
tab <- rbind(b.error, omega.error)
rownames(tab) <- c("EL Directed", "RICF Directed", "Euclidean Directed",
                   "EL Bi-directed", "RICF Bi-directed", "Euclidean Bi-directed")
colnames(tab) <- c(100, 250, 500, 1000)
xtable(tab, caption = paste("MedianSE for", type, "distributed errors in true model. All three methods converged", sum(valid), "times out of 500."), 
       digits = 4)
