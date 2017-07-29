getPolyCI <- function(y, sigma, testMat, threshold, confidence_level,
                      computeCI = TRUE) {
  etaMat <- diag(length(y))
  polyhedral <- apply(etaMat, 2, polyhedral.workhorse,
                      u = y, sigma = sigma, testMat = testMat,
                      threshold = threshold, computeCI = computeCI,
                      alpha = confidence_level)
  polyPval <- sapply(polyhedral, function(x) x$pval)
  polyCI <- do.call("rbind", lapply(polyhedral, function(x) x$ci))
  return(list(pval = polyPval, ci = polyCI))
}

getSwitchCI <- function(y, muhat, nullDist, testMat, switchTune, confidence_level, pthreshold, quadlam,
                        naiveCI) {
  p <- length(y)
  if(switchTune == "half") {
    switchLevel <- confidence_level / 2
  } else if(switchTune == "sqrd") {
    switchLevel <- confidence_level^2
  }
  t2 <- switchLevel * pthreshold
  t2threshold <- getQuadraticThreshold(t2, quadlam)
  testStat <- as.numeric(t(y) %*% testMat %*% y)
  if(testStat > t2threshold) {
    switchCI <- naiveCI
  } else {
    switchLevel <- confidence_level - switchLevel
    switchCI <- getNullCI(muhat, nullDist, switchLevel)
  }
  return(switchCI)
}

getNaiveCI <- function(y, sigma, confidence_level) {
  p <- length(y)
  sds <- sqrt(diag(sigma))
  quant <- qnorm(1 - confidence_level / 2)
  naiveCI <- matrix(nrow = p, ncol = 2)
  for(i in 1:p) {
    naiveCI[i, ] <- y[i] + c(-1, 1) * sds[i] * quant
  }
  return(naiveCI)
}

getNullCI <- function(muhat, nullDist, confidence_level) {
  p <- length(muhat)
  nullCI <- matrix(nrow = p, ncol = 2)
  for(i in 1:nrow(nullCI)) {
    nullCI[i, ] <- muhat[i] - quantile(nullDist[, i], c(1 - confidence_level / 2, confidence_level / 2))
  }

  if(nrow(nullDist) < 1 / confidence_level * 5) {
    warning("Number of samples from the truncated distribution may be insufficient
            for accurate CI computation!")
  }

  return(nullCI)
}

getMleCI <- function(muhat, mleDist, confidence_level) {
  p <- length(muhat)
  mleCI <- matrix(nrow = p, ncol = 2)
  q <- confidence_level / 2
  for(i in 1:p) {
    mleCI[i, ] <- muhat[i] - quantile(mleDist[, i], c(1 - q, q))
  }

  if(nrow(mleDist) < 1 / confidence_level * 5) {
    warning("Number of samples from the truncated distribution may be insufficient
            for accurate CI computation!")
  }
  return(mleCI)
}
