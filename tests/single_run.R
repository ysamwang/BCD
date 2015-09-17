j <- 1
p <- paramSettings[j,1]
n <-paramSettings[j,2]
k <- paramSettings[j,3]
d <- paramSettings[j,4]
b <-paramSettings[j,4]/2


out.ricf <- ricf(B = B, Omega = Omega, Y = Y, BInit = NULL,
                 OmegaInit = NULL, sigConv = 0, maxIter =5000,
                 msgs = FALSE, omegaInitScale = .9)

out.sem <- sem(model = specifyModelSEM(B, Omega), data = d)
dat <- as.data.frame(t(Y))
names(dat) <- paste("x", c(1:p), sep = "")
time.sem <- try(microbenchmark::microbenchmark(out.sem <- sem(model = specifyModelSEM(B, Omega), data = dat), times = 1),
                silent = T)