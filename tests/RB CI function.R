# Generating Data -------------
n <- 1
p <- 20
sparsity <- 2
snr <- 0
mu <- c(rep(snr / sparsity, sparsity), rep(0, p - sparsity))
sigma <- cov2cor(clusterGeneration::genPositiveDefMat(p)$Sigma)
sqrtSig <- expm::sqrtm(sigma)
testmat <- solve(sigma) * n

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
y <- as.numeric(y)

# Computing matrices ---------------
testEigen <- eigen(testmat)
invTestMat <- testEigen$vectors %*% diag(1/(testEigen$values))
invTestMat <- invTestMat %*% t(testEigen$vectors)
sqrtTestMat <- testEigen$vectors %*% diag(sqrt(testEigen$values))
sqrtTestMat <- sqrtTestMat %*% t(testEigen$vectors)
invSqrtTestMat <- testEigen$vectors %*% diag(1 / sqrt(testEigen$values))
invSqrtTestMat <- invSqrtTestMat %*% t(testEigen$vectors)
sampSig <- t(sqrtTestMat) %*% sigma %*% sqrtTestMat
sampP <- solve(sampSig)

# Computing CIs -------------
alpha <- 0.05
rbIters <- (40 * (1 - alpha) / alpha)
ci <- matrix(nrow = p, ncol = 2)
for(j in 1:p) {
  for(upper in c(TRUE, FALSE)) {
    if(upper) {
      init <- y[j] + qnorm(1 - alpha / 4) * sqrt(diag(sigma)[j])
      a <- alpha / 2
    } else {
      init <- y[j] - qnorm(1 - alpha / 4) * sqrt(diag(sigma)[j])
      a <- alpha / 2
    }

    u <- quadraticRobinsMonroe(var = j, upper = upper, alpha = a,
                               observed = y[j], initU = init,
                               threshold = threshold,
                               sqrtTestMat = sqrtTestMat,
                               invSqrtTestMat = invSqrtTestMat,
                               sampPrecision = sampP,
                               stepSize = 17 * sqrt(diag(sigma)[j]),
                               burnin = 1000, mhiters = 5, rbIters = rbIters)
    ci[j, upper + 1] <- u
  }
  print(ci[j, ])
}


maxtries <- 50
uvec <- numeric(maxtries)
for(try in 1:maxtries) {
  uvec[try] <- u
  print(u)
  if(try > 1) {
    ratio <- abs(1 - u/oldu)
    # print(ratio)
    # if(ratio < 10^-2) break
  }
  oldu <- u
}
u <- (u + oldu) / 2
print(mean(uvec))
print(sd(uvec))
hist(uvec)
# par(mfrow = c(3, 1), mar = rep(3, 4))




