library(ggplot2)
library(dplyr)
library(reshape2)
n <- 1
p <- 20
sparsity <- 2
snr <- 0
mu <- c(rep(snr / sparsity, sparsity), rep(0, p - sparsity))
sigma <- cov2cor(clusterGeneration::genPositiveDefMat(p)$Sigma)
sqrtSig <- expm::sqrtm(sigma)
testmat <- solve(sigma) * n

# sigma <- diag(p)
pthreshold <- 0.05
threshold <- qchisq(pthreshold, p, lower.tail = FALSE)
testStat <- 0
invcov <- solve(sigma)
while(testStat < threshold) {
  z <- matrix(rnorm(n * p), nrow = n)
  y <- z %*% sqrtSig
  naive <- colMeans(y)
  testStat <- as.numeric(t(naive) %*% testmat %*% naive)
}

# Testing line-search function --------
ncp <- as.numeric(naive %*% testmat %*% naive)
optimLambda <- optimize(conditionalDnorm, interval = c(0, 1),
                        maximum = TRUE, y = naive, precision = invcov * n,
                        ncp = ncp, threshold = threshold)$maximum
muhat <- optimLambda * naive

# Testing SGD function -------------
nsteps <- 1000
sgdfit <- quadraticSGD(naive, sigma / n, invcov * n, testmat, threshold,
                      stepRate = 0.85, stepCoef = NULL,
                      delay = 10, sgdSteps = nsteps, assumeCovergence = nsteps / 2,
                      mhIters = 40)
solutionPath <- sgdfit$solutionPath
muSGD <- solutionPath[nrow(solutionPath), ]
cor(muhat, muSGD)
plot(naive, muhat, col = "blue", ylim = c(min(muhat, muSGD), max(muhat, muSGD)))
points(naive, muSGD, col = "red")
abline(a = 0, b = 1)

solutionPath <- melt(solutionPath)
names(solutionPath) <- c("iter", "var", "muhat")
ggplot(solutionPath) + geom_line(aes(x = iter, y = muhat, col = factor(var))) +
  theme_bw() + geom_hline(yintercept = 0)

# Testing null sampler
confidence_level <- 0.05
nSamples <- 10 / confidence_level
nSamples <- max(nSamples, 30 * p, 1000)
nullMu <- rep(0, p)
nullSample <- sampleQuadraticConstraint(nullMu, sigma / n,
                                        threshold, testmat,
                                        sampSize = nSamples,
                                        burnin = 1000,
                                        trim = 50)
nullPval <- numeric(p)
for(i in 1:length(nullPval)) {
  nullPval[i] <- 2 * min(mean(naive[i] < nullSample[, i]), mean(naive[i] > nullSample[, i]))
}
naivePval <- 2*pnorm(-abs(naive), sd = sqrt(diag(sigma) / n))
cbind(naivePval, nullPval)

# Testing polyhedral functions ----------------------
etaMat <- diag(p)
computeCI <- TRUE
polyhedral <- apply(etaMat, 2, polyhedral.workhorse,
                    u = naive, sigma = sigma / n, testMat = testmat,
                    threshold = threshold, computeCI = computeCI,
                    alpha = 0.05)
polyPval <- sapply(polyhedral, function(x) x$pval)
polyCI <- do.call("rbind", lapply(polyhedral, function(x) x$ci))
hybridPval <- pmin(2 * pmin(polyPval, nullPval), 1)
cbind(naivePval, nullPval, polyPval, hybridPval)

# MLE based inference ----------
mleSample <- sampleQuadraticConstraint(muhat, sigma / n,
                                       threshold, testmat,
                                       sampSize = nSamples,
                                       burnin = 1000,
                                       trim = 100)
if(p < 20) {
  mleLam <- apply(mleSample, 1, function(x) optimize(conditionalDnorm, interval = c(0, 1),
                                                     maximum = TRUE, y = x, precision = invcov,
                                                     ncp = ncp, threshold = threshold)$maximum)
} else {
  mleLam <- rep(1, nrow(mleSample))
}

mleDist <- mleSample
for(i in 1:nrow(mleDist)) {
  mleDist[i, ] <- mleDist[i, ] * mleLam[i] - muhat
}

mleDist <- apply(ml)
mlePval <- numeric(p)
for(i in 1:length(mlePval)) {
  if(muhat[i] > 0) {
    if(muhat[i] - max(mleDist[, i]) > 0) {
      mlePval[i] <- 1 / nrow(mleDist)
    } else {
      mlePval[i] <- 1 - uniroot(function(x) muhat[i] - quantile(mleDist[, i], x),
                                interval = c(0, 1))$root
    }
  } else {
    if(muhat[i] - min(mleDist[, i]) < 0) {
      mlePval[i] <- 1 / nrow(mleDist)
    } else {
      mlePval[i] <- uniroot(function(x) muhat[i] - quantile(mleDist[, i], x),
                            interval = c(0, 0.5))$root
    }
  }
}
mlePval <- 2 * mlePval
cbind(naivePval, nullPval, polyPval, hybridPval, mlePval)

mleCI <- matrix(nrow = p, ncol = 2)
for(i in 1:nrow(mleCI)) {
  mleCI[i, ] <- muhat[i] - quantile(mleDist[, i], c(1 - confidence_level / 2, confidence_level / 2))
}
cbind(polyCI, mleCI)

# Global null CIs ---------------
if(p < 20) {
  nullLam <- apply(nullSample, 1, function(x) optimize(conditionalDnorm, interval = c(0, 1),
                                                       maximum = TRUE, y = x, precision = invcov * n,
                                                       ncp = ncp, threshold = threshold)$maximum)
} else {
  nullLam <- rep(optimLambda, nrow(nullSample))
}

nullDist <- nullSample
for(i in 1:nrow(nullDist)) {
  nullDist[i, ] <- nullDist[i, ] * nullLam[i]
}

nullCI <- matrix(nrow = p, ncol = 2)
for(i in 1:nrow(nullCI)) {
  nullCI[i, ] <- muhat[i] - quantile(nullDist[, i], c(1 - confidence_level / 2, confidence_level / 2))
}

nullCI
