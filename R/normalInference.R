mvnQuadratic <- function(y, sigma, testMat = "wald",
                         threshold = NULL, pval_threshold = 0.05,
                         estimate_type = c("mle", "naive"),
                         pvalue_type = c("polyhedral", "hybrid", "naive", "global-null", "mle"),
                         ci_type = c("polyhedral", "switch", "naive", "global-null", "mle"),
                         confidence_level = 0.05,
                         switchTune = c("sqrd", "half"),
                         nSamples = NULL, verbose = TRUE) {
  # Validating parameters --------------------
  if(!(any(estimate_type %in% c("mle", "naive")))) {
    stop("estimation method not supported!")
  }

  if(!(any(pvalue_type %in% c("polyhedral", "hybrid", "global-null", "mle")))) {
    stop("pvalue_type not supported!")
  }

  if(!(any(ci_type %in% c("polyhedral", "switch", "naive", "mle")))) {
    stop("ci_type not supported!")
  }

  if("mle" %in% ci_type & !("mle" %in% estimate_type)) {
    warning("Computing the conditional MLE is necessary for computing MLE based confidence intervals.")
    estimate_type <- c(estimate_type, "mle")
  }

  if(("mle" %in% c(ci_type, pvalue_type)) & length(y) >= 10) {
    warning("mle based inference is not recommended for dim >= 10!")
    if(testMat[1] != "wald") {
      warning("MLE based inference is conservative for non-Wald tests!")
    }
  }

  if(testMat[1] == "wald") {
    quadtest <- "wald"
    invcov <- solve(sigma)
    testMat <- invcov
    quadlam <- rep(1, p)
  } else {
    invcov <- solve(sigma)
    quadtest <- "other"
    quadlam <- getQudraticLam(testMat, sigma)
  }

  if(any(length(y) != dim(sigma))) {
    stop("Incorrect dimensions for y or sigma!")
  }

  if(any(dim(testMat) != dim(sigma))) {
    stop("Dimension of sigma must match the dimesions of testMat!")
  }

  if(confidence_level < 0 | confidence_level > 1) {
    stop("confidence_level must be between 0 and 1!")
  }

  switchTune <- switchTune[1]
  if(!(switchTune %in% c("sqrd", "half"))) {
    stop("Regime switching CI tuning parameter type not supported!")
  }

  # Setting threshold ----------------------------
  p <- nrow(sigma)
  if(is.null(threshold)) {
    if(pthreshold < 0 | pthreshold > 1) {
      stop("If a threshold is not provided, then pthreshold must be between 0 and 1!")
    } else {
      threshold <- getQuadraticThreshold(pthreshold, quadlam)
    }
  }

  pthreshold <- CompQuadForm::liu(threshold, quadlam)

  # Computing the MLE ---------------------
  if(verbose) {
    print("Computing MLE!")
  }

  if(estimate_type[1] == "naive") {
    muhat <- y
    solutionPath <- NULL
    mleMu <- NULL
  }

  if( "mle" %in% estimate_type) {
    if(quadtest == "wald") {
      ncp <- as.numeric(t(y) %*% invcov %*% y)
      optimLambda <- optimize(conditionalDnorm, interval = c(0, 1),
                              maximum = TRUE, y = y, precision = invcov,
                              ncp = ncp, threshold = threshold)$maximum
      mleMu <- optimLambda * y
      solutionPath <- NULL
    } else if(quadtest == "other") {
      sgdfit <- quadraticSGD(y, sigma, invcov, testMat, threshold,
                             stepRate = 0.75, stepCoef = NULL,
                             delay = 20, sgdSteps = 1000, assumeCovergence = 800,
                             mhIters = 40)
      mleMu <- sgdfit$mle
      solutionPath <- sgdfit$solutionPath
      optimLambda <- NULL
    }
  } else {
    mleMu <- NULL
    solutionPath <- NULL
    optimLambda <- NULL
  }

  if(estimate_type[1] == "mle") {
    muhat <- mleMu
  }

  # Sampling from the global-null ---------------------
  if(is.null(nSamples)) {
    nSamples <- length(y) * 5 / confidence_level
    nSamples <- max(nSamples, 30 * length(y))
  }

  if("global-null" %in% pvalue_type | "switch" %in% ci_type) {
    if(verbose) print(paste("Sampling from null distribution! (", nSamples, " samples)", sep = ""))
    nullMu <- rep(0, length(y))
    nullSample <- sampleQuadraticConstraint(nullMu, sigma,
                                            threshold, testMat,
                                            sampSize = nSamples,
                                            burnin = 1000,
                                            trim = 50)
    nullPval <- numeric(length(y))
    for(i in 1:length(nullPval)) {
      nullPval[i] <- 2 * min(mean(y[i] < nullSample[, i]), mean(y[i] > nullSample[, i]))
    }
    if(pvalue_type[1] == "global-null") {
      pvalue <- nullPval
    }
  } else {
    nullSample <- NULL
    nullPval <- NULL
    nullDist <- NULL
  }

  # Computing polyhedral p-values/CIs ---------------------
  if(any(c("polyhedral", "hybrid") %in% pvalue_type) | "polyhedral" %in% ci_type) {
    if(verbose) print("Computing polyhedral p-values/CIs!")
    polyResult <- getPolyCI(y, sigma, testMat, threshold, confidence_level)
    polyPval <- polyResult$pval
    polyCI <- polyResult$ci
    if(pvalue_type[1] == "polyhedral") pvalue <- polyPval
    if(ci_type[1] == "polyhedral") ci <- polyCI
  } else {
    polyPval <- NULL
    polyCI <- NULL
  }

  # Computing hybrid p-value -----------------
  if("hybrid" %in% pvalue_type) {
    hybridPval <- pmin(2 * pmin(polyPval, nullPval), 1)
  } else {
    hybridPval <- NULL
  }

  if(pvalue_type[1] == "hybrid") {
    pvalue <- hybridPval
  }

  # MLE based inference -----------
  if("mle" %in% c(pvalue_type, ci_type)) {
    if(verbose) print(paste("Sampling from the truncated distribution under the MLE! (", nSamples, " samples)", sep =""))
    mleSample <- sampleQuadraticConstraint(muhat, sigma,
                                           threshold, testMat,
                                           sampSize = nSamples,
                                           burnin = 1000,
                                           trim = 50)
    if(p < 20 & quadtest == "wald") {
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
                                interval = c(0, 1))$root
        }
      }
    }
    mlePval <- 2 * mlePval
    mleCI <- getMleCI(muhat, mleDist, confidence_level)
  } else {
    mleSample <- NULL
    mlePval <- NULL
    mleCI <- NULL
    mleDist <- NULL
  }

  if(ci_type[1] == "mle") {
    ci <- mleCI
  }
  if(pvalue_type[1] == "mle") {
    pvalue <- mlePval
  }

  # Null Distribution Based CIs ---------------
  if(any(c("switch", "global-null") %in% ci_type)) {
    if(p < 20 & quadtest == "wald") {
      if(verbose) print("Bootstrapping the null-distribution of the MLE under the null!")
      nullLam <- apply(nullSample, 1, function(x) optimize(conditionalDnorm, interval = c(0, 1),
                                                           maximum = TRUE, y = x, precision = invcov,
                                                           ncp = ncp, threshold = threshold)$maximum)
    } else {
      if(p < 20 & quadtest == "other") {
        warning("Switching regime and null distirbution based CIs may be conservative for general qudaratic tests!")
      }
      nullLam <- rep(1, nrow(nullSample))
    }

    nullDist <- nullSample
    for(i in 1:nrow(nullDist)) {
      nullDist[i, ] <- nullDist[i, ] * nullLam[i]
    }
    nullCI <- getNullCI(muhat, nullDist, confidence_level)
  } else {
    nullDist <- NULL
    nullCI <- NULL
  }

  if(ci_type[1] == "global_null") {
    ci <- nullCI
  }

  # Naive pvalues and CIs -----------------------
  naiveCI <- getNaiveCI(y, sigma, confidence_level)
  naivePval <- 2 * pnorm(-abs(as.numeric(y / sqrt(diag(sigma)))))
  if(ci_type[1] == "naive") {
    ci <- naiveCI
  }
  if(pvalue_type[1] == "naive") {
    pvalue <- naivePval
  }

  # Regime switching CIs -------------------------
  if(any(c("switch", "global-null") %in% ci_type)) {
    switchCI <- getSwitchCI(y, muhat, nullDist, testMat, switchTune,
                            confidence_level, pthreshold, quadlam,
                            naiveCI)
  } else {
    switchCI <- NULL
  }

  if(ci_type[1] == "switch") {
    ci <- switchCI
  }

  # Wrapping up ------------------
  results <- list()
  results$muhat <- muhat
  results$ci <- ci
  results$pvalue <- pvalue

  results$naiveMu <- y
  results$mleMu <- mleMu
  results$solutionPath <- solutionPath
  results$nullSample <- nullSample
  results$nullDist <- nullDist
  results$nullPval <- nullPval
  results$polyPval <- polyPval
  results$polyCI <- polyCI
  results$mleSample <- mleSample
  results$mleDist <- mleDist
  results$mlePval <- mlePval
  results$mleCI <- mleCI
  results$nullDist <- nullDist
  results$nullCI <- nullCI
  results$naivePval <- naivePval
  results$naiveCI <- naiveCI
  results$switchCI <- switchCI
  results$hybridPval <- hybridPval

  results$testMat <- testMat
  results$sigma <- sigma
  results$invcov <- invcov
  results$pthreshold <- pthreshold
  results$threshold <- threshold
  results$confidence_level <- confidence_level
  results$switchTune <- switchTune
  results$estimate_type <- estimate_type
  results$pvalue_type <- pvalue_type
  results$ci_type <- ci_type

  results$quadlam <- quadlam

  class(results) <- "mvnQuadratic"
  if(verbose) print("Done!")
  return(results)
}

conditionalDnorm <- function(lambda, y, precision, ncp, threshold) {
  mu <- y * lambda
  ncp <- ncp * lambda^2

  prob <- pchisq(threshold, length(y), ncp, FALSE, TRUE)

  diff <- y - mu
  density <- -0.5 * t(diff) %*% precision %*% diff
  condDens <- density - prob
  return(condDens)
}

getQudraticLam <- function(testMat, sigma) {
  c <- chol(sigma)
  tempmat <- c %*% testMat %*% t(c)
  eig <- eigen(tempmat)
  vec <- eig$vectors
  P <- t(vec)
  lam <- eig$values
  return(lam)
}

getQuadraticThreshold <- function(quantile, lam) {
  lower <- 0
  upper <- 1
  prob <- CompQuadForm::liu(upper, lam)
  while(prob > quantile) {
    lower <- upper
    upper <- upper * 2
    prob <- CompQuadForm::liu(upper, lam)
  }

  q <- uniroot(function(x) CompQuadForm::liu(x, lam) - quantile,
               interval = c(lower, upper))$root
  return(q)
}
