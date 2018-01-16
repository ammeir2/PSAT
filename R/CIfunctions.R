getPolyCI <- function(y, sigma, testMat, threshold, confidence_level,
                      computeCI = TRUE, test = "quadratic",
                      truncPmethod = "symmetric") {
  etaMat <- diag(length(y))
  if(test == "quadratic") { # Quadratic
    polyhedral <- apply(etaMat, 2, polyhedral.workhorse,
                        u = y, sigma = sigma, testMat = testMat,
                        threshold = threshold, computeCI = computeCI,
                        alpha = confidence_level, truncPmethod = truncPmethod)
    noCorrection <- sapply(polyhedral, function(x) x$noCorrection)
  } else if(test == "linear") { # Linear
    polyhedral <- apply(etaMat, 2, linearPolyhedral,
                        u = y, sigma = sigma, a = testMat,
                        threshold = threshold, computeCI = computeCI,
                        alpha = confidence_level, truncPmethod = truncPmethod)
    noCorrection <- rep(FALSE, length(y))
  }
  polyPval <- sapply(polyhedral, function(x) x$pval)
  polyCI <- do.call("rbind", lapply(polyhedral, function(x) x$ci))
  return(list(pval = polyPval, ci = polyCI, noCorrection = noCorrection))
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

getHybridCI <- function(y, sigma, testMat, threshold, pthreshold, confidence_level,
                        hybridPval = NULL, trueHybrid = FALSE, rbIters = NULL,
                        test = "quadratic") {
  ci <- getSwitchCI(y, sigma, testMat, threshold, pthreshold, confidence_level,
                    quadlam = rep(1, length(y)),
                    t2 = 10^-10,
                    testStat = 0, rbIters = rbIters, test)
  return(ci)
}

getSwitchCI <- function(y, sigma, testMat, threshold, pthreshold,
                        confidence_level,
                        quadlam, t2, testStat, hybridPval = NULL,
                        trueHybrid = FALSE, rbIters = NULL,
                        test = "quadratic") {
  if(test == "linear") {
    contrast <- testMat
  }

  p <- length(y)
  # Place holder for hybridPval
  if(is.null(hybridPval)) {
    hybridPval <- rep(0, p)
  }

  # Naive or Hybrid?
  if(test == "quadratic") {
    secondThreshold <- getQuadraticThreshold(t2, quadlam)
    if(testStat > secondThreshold) {
      ci <- t(sapply(1:length(y), function(j) y[j] + c(-1, 1) * qnorm(1 - confidence_level / 2) * sqrt(sigma[j, j])))
      return(ci)
    }
  } else if(test == "linear") {
    contsd <- as.numeric(sqrt(t(contrast) %*% sigma %*% contrast))
    secondThreshold <- numeric(2)
    secondThreshold[1] <- qnorm(pnorm(threshold[1], sd = contsd) * t2, sd = contsd)
    secondThreshold[2] <- qnorm(pnorm(threshold[2], sd = contsd, lower.tail = FALSE) * t2,
                                sd = contsd, lower.tail = FALSE)
    testStat <- sum(y * contrast)
    if(testStat < secondThreshold[1] | testStat > secondThreshold[2]) {
      ci <- t(sapply(1:length(y), function(j) y[j] + c(-1, 1) * qnorm(1 - confidence_level / 2) * sqrt(sigma[j, j])))
      return(ci)
    }
  }

  # Computing polyhedral CIs ----------
  if(test == "quadratic") {
    a <- (confidence_level - t2 * pthreshold) / 2
  } else if(test == "linear") {
    a <- (confidence_level - t2) / 2
  }
  polyCI <- getPolyCI(y, sigma, testMat, threshold, a , computeCI = TRUE,
                      test = test)
  poly <- polyCI$ci

  if(!trueHybrid) {
    noBonferroni <- getPolyCI(y, sigma, testMat, threshold, 2 * a, computeCI = TRUE,
                              test = test)$ci
    poly[y < 0, 1] <- noBonferroni[y < 0, 1]
    poly[y > 0, 2] <- noBonferroni[y > 0, 2]
  }

  # Computing global-null CIs ---------
  if(trueHybrid) {
    whichCompute <- which(!polyCI$noCorrection)
    computeFull <- TRUE
  } else {
    whichPval <- (hybridPval < confidence_level) & (polyCI$pval > confidence_level)
    whichCompute <- which(!polyCI$noCorrection & whichPval)
    computeFull <- FALSE
  }

  if(length(whichCompute) > 0) {
    if(test == "quadratic") {
      nullCI <- quadraticRB(y, sigma, testMat, threshold, a,
                            variables = whichCompute,
                            computeFull = computeFull,
                            rbIters = rbIters)
    } else if(test == "linear") {
      nullCI <- linearRB(y, sigma, contrast, threshold, a,
                         rbIters = rbIters,
                         variables = whichCompute,
                         computeFull = computeFull)
    }
  } else {
    nullCI <- cbind(rep(-Inf, p), rep(Inf, p))
  }

  ci <- cbind(pmax(poly[, 1], nullCI[, 1]),
              pmin(poly[, 2], nullCI[, 2]))

  return(ci)
}
