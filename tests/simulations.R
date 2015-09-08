#### Simulation Studies for Paper ####
library(MASS)
library(BCD)
library(sem)
library(microbenchmark)
source("tests/simHelper.R")

set.seed(1001)
p <- 10
n <- 500
k <- 5
d <- .2
b <- d/2

ret <- do.one(p, n, k, d, b, times = 10)
ret$ricfTime/ret$semTime
ret$agree
ret$ricfConv
ret$semConv