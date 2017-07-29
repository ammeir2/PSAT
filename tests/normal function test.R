# Tests for wald test ------------------------
n <- 1
p <- 6
sparsity <- 2
snr <- 3
mu <- c(rep(snr / sparsity, sparsity), rep(0, p - sparsity))
sigma <- cov2cor(clusterGeneration::genPositiveDefMat(p)$Sigma)
sqrtSig <- expm::sqrtm(sigma)
testmat <- solve(sigma) * n

pthreshold <- 0.001
threshold <- qchisq(pthreshold, p, lower.tail = FALSE)
testStat <- 0
invcov <- solve(sigma)
while(testStat < threshold) {
  z <- matrix(rnorm(n * p), nrow = n)
  y <- as.numeric(z %*% sqrtSig) + mu
  testStat <- as.numeric(t(y) %*% testmat %*% y)
}

allfit <- mvnQuadratic(y, sigma, testMat = testmat, threshold = NULL, pval_threshold = pthreshold,
                       estimate_type = c("mle", "naive"),
                       pvalue_type = c("polyhedral", "hybrid", "naive", "global-null"),
                       ci_type = c("polyhedral", "switch", "naive", "global-null"),
                       confidence_level = 0.05,
                       switchTune = c("sqrd", "half"),
                       nSamples = NULL, verbose = TRUE)

plot(allfit, type = "estimates", true = mu)

fitlist <- list()
slot <- 1
pvaltypes <- c("polyhedral", "hybrid", "global-null", "mle")
citypes <- c("polyhedral", "switch", "naive", "mle")
esttypes <- c("mle", "naive")
ntests <- length(pvaltypes) * length(citypes) * length(esttypes)
for(i in 1:length(pvaltypes)) {
  for(j in 1:length(citypes)) {
    for(k in 1:length(esttypes)) {
      fitlist[[slot]] <- mvnQuadratic(y, sigma, testMat = "wald",
                                      pval_threshold = 0.05, confidence_level = 0.05,
                                      estimate_type = esttypes[k],
                                      pvalue_type = pvaltypes[i],
                                      ci_type = citypes[j],
                                      switchTune = c("sqrd"),
                                      nSamples = NULL, verbose = FALSE)
      cat(round(slot / ntests, 2)*100, "% ")
      slot <- slot + 1
    }
  }
}

# Tests for general quadratic test ------------------------
n <- 1
p <- 5
sparsity <- 2
snr <- 0
mu <- c(rep(snr / sparsity, sparsity), rep(0, p - sparsity))
sigma <- cov2cor(clusterGeneration::genPositiveDefMat(p)$Sigma)
sqrtSig <- expm::sqrtm(sigma)
testmat <- cov2cor(clusterGeneration::genPositiveDefMat(p)$Sigma)
# testmat <- solve(sigma)
quadlam <- getQudraticLam(testmat, sigma)
pthreshold <- 0.05
threshold <- getQuadraticThreshold(pthreshold, quadlam)
testStat <- 0
invcov <- solve(sigma)
while(testStat < threshold) {
  z <- matrix(rnorm(n * p), nrow = n)
  y <- as.numeric(z %*% sqrtSig)
  testStat <- as.numeric(t(y) %*% testmat %*% y)
}

fitlist <- list()
slot <- 1
pvaltypes <- c("polyhedral", "hybrid", "global-null", "mle")
citypes <- c("polyhedral", "switch", "naive", "mle")
esttypes <- c("mle", "naive")
ntests <- length(pvaltypes) * length(citypes) * length(esttypes)
for(i in 1:length(pvaltypes)) {
  for(j in 1:length(citypes)) {
    for(k in 1:length(esttypes)) {
      fitlist[[slot]] <- mvnQuadratic(y, sigma, testMat = testmat,
                                      pval_threshold = 0.05, confidence_level = 0.05,
                                      estimate_type = esttypes[k],
                                      pvalue_type = pvaltypes[i],
                                      ci_type = citypes[j],
                                      switchTune = c("sqrd"),
                                      nSamples = NULL, verbose = FALSE)
      cat(round(slot / ntests, 2)*100, "% ")
      slot <- slot + 1
    }
  }
}

fitlist[[1]]$muhat / fitlist[[1]]$naiveMu
