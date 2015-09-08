#### Simulation Studies for Paper ####
library(MASS)
library(BCD)
library(sem)
library(microbenchmark)
source("simHelper.R")


v <- c(rep(10, 12), rep(20, 12))
n <- c(rep(15, 6), rep(100, 6), rep(30, 6), rep(200, 6))
k <- c(rep(rep(c(0,2,4), each = 2),2), rep(rep(c(0,4,8), each = 2),2))
d <- rep(c(.1, .2), 12)
paramSettings <- data.frame(v = v, n = n, k = k, d = d)    

set.seed(1001)


sim.size <- 1000
ricfTime <- semTime <- ricfConv <- semConv <- agree <- ricfIter <- ricfCrit <- rep(0, sim.size)
j <- 1

for(i in 1:sim.size){
  ret <- do.one(paramSettings[j,1], paramSettings[j,2], paramSettings[j,3], paramSettings[j,4], paramSettings[j,4]/2, times = 1)
  ricfTime[i] <- ret$ricfTime
  semTime[i] <- ret$semTime
  ricfConv[i] <- ret$ricfConv
  semConv[i] <- ret$semConv
  agree[i] <- ret$agree
}

mean(ricfConv)
mean(semConv)
mean(agree[agree >0])
