#### Simulation Studies for Paper ####
source("tests/simHelper.R")


v <- c(rep(10, 12), rep(20, 12))
n <- c(rep(15, 6), rep(100, 6), rep(30, 6), rep(200, 6))
k <- c(rep(rep(c(0,2,4), each = 2),2), rep(rep(c(0,4,8), each = 2),2))
d <- rep(c(.1, .2), 12)
paramSettings <- data.frame(v = v, n = n, k = k, d = d)    
results <- data.frame(bcdConv = rep(0, dim(paramSettings)[1]),
                      semConv = rep(0, dim(paramSettings)[1]),
                      mixedConv = rep(0, dim(paramSettings)[1]),
                      bothConv = rep(0, dim(paramSettings)[1]),
                      agree = rep(0, dim(paramSettings)[1]),
                      bcdTime = rep(0, dim(paramSettings)[1]),
                      semTime = rep(0, dim(paramSettings)[1]),
                      mixedTime = rep(0, dim(paramSettings)[1]),
                      bcdTimeComparison = rep(0, dim(paramSettings)[1]),
                      semTimeComparison = rep(0, dim(paramSettings)[1]),
                      mixedTimeComparison = rep(0, dim(paramSettings)[1]))

set.seed(111)


sim.size <- 1000


for(j in 1:dim(paramSettings)[1]){
  print(paste("Setting", j, sep = " "))
  print(paramSettings[j,])
  
  ricfTime <- semTime <- mixedTime <- ricfConv <- semConv <- mixedConv <- agree <- rep(0, sim.size)
  for(i in 1:sim.size){
    cat(paste (".", ifelse(i%%50==0,i,""), sep = ""))
    ret <- do.one(paramSettings[j,1], paramSettings[j,2], paramSettings[j,3], paramSettings[j,4], paramSettings[j,4]/2, times = 1)
    ricfTime[i] <- ret$ricfTime
    semTime[i] <- ret$semTime
    mixedTime[i] <- ret$mixedTime
    ricfConv[i] <- ret$ricfConv
    semConv[i] <- ret$semConv
    mixedConv[i] <- ret$mixedConv
    agree[i] <- ret$agree
  }

  results[j, ] <- c(sum(ricfConv),
                    sum(semConv),
                    sum(mixedConv),
                    sum(ricfConv * semConv),
                    sum(agree),
                    mean(ricfTime[ricfConv == 1]),
                    mean(semTime[semConv == 1]),
                    mean(mixedTime[mixedConv == 1]),
                    mean(ricfTime[agree == 1]),
                    mean(semTime[agree == 1]),
                    mean(mixedTime[agree == 1]))
  cat("\n")
}

# saveRDS(results, "data/table3.RDS")

results <- readRDS("data/table3.RDS")
results[,6:11] <- results[,6:11] * 1e-6
tab <- cbind(paramSettings, results[-c(3,6:8,11)])
library(xtable)
print(xtable(tab, digits = c(0,0,0,0,1,0,0,0,0,1,1)), include.rownames =F)
