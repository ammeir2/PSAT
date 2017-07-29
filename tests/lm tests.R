n <- 500
p <- 50
sparsity <- 3
snr <- 0.1
beta <- rnorm(sparsity)
beta <- beta / sum(abs(beta)) * snr
beta <- c(beta, rep(0, p - sparsity))
xcov <- cov2cor(clusterGeneration::genPositiveDefMat(p)$Sigma)
pthreshold <- 0.05
threshold <- qchisq(1 - pthreshold, p)
X <- as.matrix(scale(mvtnorm::rmvnorm(n, sigma = xcov)))
mu <- as.numeric(X %*% beta)
testMat <- t(X) %*% X
projmat <- solve(t(X) %*% X) %*% t(X)
testStat <- 0
while(testStat < threshold) {
  y <- rnorm(n) + mu
  y <- y - mean(y)
  naive <- as.numeric(projmat %*% y)
  testStat <- as.numeric(t(naive) %*% testMat %*% naive) / sd(y)
}

fit <- lmQuadratic(X, y, testMat = "wald", resid_sd = "ysd",
                   threshold = NULL, pval_threshold = pthreshold,
                   estimate_type = "mle",
                   pvalue_type = "hybrid",
                   ci_type = "polyhedral",
                   confidence_level = 0.05,
                   switchTune = c("sqrd", "half"),
                   nSamples = NULL, verbose = TRUE)


