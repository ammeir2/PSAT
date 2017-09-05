quadraticSGD <- function(y, sigma, precision, testmat, threshold,
                         stepRate = 0.75, stepCoef = NULL,
                         delay = 10,
                         sgdSteps = 800, assumeCovergence = 400,
                         mhIters = 40) {

  if(is.null(stepCoef)) {
    stepCoef <- 0.25 / sqrt(diag(sigma))
  } else if(length(stepCoef) == 1) {
    stepCoef <- stepCoef / sqrt(diag(sigma))
  }

  # Pre-processing test matrix ---------
  testEigen <- eigen(testmat)
  invTestMat <- testEigen$vectors %*% diag(1/(testEigen$values))
  invTestMat <- invTestMat %*% t(testEigen$vectors)
  sqrtTestMat <- testEigen$vectors %*% diag(sqrt(testEigen$values))
  sqrtTestMat <- sqrtTestMat %*% t(testEigen$vectors)
  invSqrtTestMat <- testEigen$vectors %*% diag(1 / sqrt(testEigen$values))
  invSqrtTestMat <- invSqrtTestMat %*% t(testEigen$vectors)
  sampSig <- t(sqrtTestMat) %*% sigma %*% sqrtTestMat
  sampP <- solve(sampSig)

  # Setting up SGD ----------------
  itermu <- y
  mu <- y
  mu <- sqrtTestMat %*% mu
  lastsamp <- as.numeric(sqrtTestMat %*% y)
  solutionPath <- matrix(ncol = length(mu), nrow = sgdSteps)
  sampmat <- matrix(ncol = length(y), nrow = sgdSteps)
  for(i in 1:sgdSteps) {
    # Sampling -----------------
    sampMu <- sqrtTestMat %*% itermu
    sample <- sampleQuadraticMVTcpp(lastsamp, sampMu, sampP, threshold,
                                    2, round(mhIters) / 2, round(mhIters) / 2)
    lastsamp <- sample[2, ]

    # computing gradient ------------
    grad <- sample %*% invSqrtTestMat
    grad <- colMeans(grad)
    grad <- precision %*% (y - grad) * stepCoef / max(1, i - delay)^stepRate
    grad <- sign(grad) * pmin(abs(grad), sqrt(diag(sigma)) * 0.05)

    # Updating estimate --------------
    itermu <- itermu + grad
    itermu <- pmax(0, itermu * sign(y)) * sign(y)
    itermu <- pmin(abs(y), abs(itermu)) * sign(y)
    solutionPath[i, ] <- itermu
    sampmat[i, ] <- lastsamp
  }

  # Reporting ----------------------
  mle <- colMeans(solutionPath[assumeCovergence:sgdSteps, ])
  return(list(mle = mle, solutionPath = solutionPath, sampMat = sampmat))
}
