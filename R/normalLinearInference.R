mvnLinear <- function(y, sigma, contrast,
                      threshold, pval_threshold,
                      selection = c("lower", "upper", "two_sided"),
                      estimate_type = c("mle", "naive"),
                      pvalue_type = c("hybrid", "polyhedral", "naive"),
                      ci_type = c("switch", "polyhedral", "naive"),
                      confidence_level = .95,
                      verbose = TRUE,
                      nSamples = NULL, trueHybrid = FALSE) {


  checkPvalues(pvalue_type)
  checkCI(ci_type)

  if(any(length(y) != dim(sigma))) {
    stop("Incorrect dimensions for y or sigma!")
  }

  if(length(y) != length(contrast)) {
    stop("Incorrect dimension for y or contrast!")
  }

  if(confidence_level <= 0 | confidence_level >= 1) {
    stop("confidence_level must be between 0 and 1!")
  }
  confidence_level <- 1 - confidence_level

  # Setting threshold ----------------------------
  if(!any(selection %in% c("two-sided", "lower", "upper")) &
     length(threshold) < 2) {
    stop("Selection must be one of `two-sided', `lower' or `upper'!
         (or a threshold of length 2 must be provided)")
  }

  p <- length(y)
  contrastVar <- as.numeric(t(contrast) %*% sigma %*% contrast)
  selection <- selection[1]
  if(is.null(threshold)) {
    if(pthreshold < 0 | pthreshold > 1) {
      stop("If a threshold is not provided, then pval_threshold must be between 0 and 1!")
    } else {
      if(selection == "two-sided") {
        threshold <- qnorm(1 - pval_threshold/2, sd = sqrt(contrastVar))
        threshold <- c(-threshold, threshold)
      } else if(selection == "lower") {
        threshold <- qnorm(pval_threshold, sd = sqrt(contrastVar))
        threshold <- c(threshold, Inf)
      } else if(selection == "upper") {
        threshold <- qnorm(1 - pval_threshold, sd = sqrt(contrastVar))
        threshold <- c(-Inf, threshold)
      } else {
        stop("If threshold is not provided then selection must be one of `two-sided'
             `lower' or `upper'!")
      }
    }
  }

  if(length(threshold) > 2) {
    stop("The length of threshold must be either 1 or 2!")
  } else if(length(threshold) == 1) {
    if(selection == "two-sided") {
      warning(paste("Threshold of length one provided for a two sided selection rule.
                    Selection rule assumed to be contrast < -abs(threshold) or contrast > abs(threshold)"))
      threshold <- c(-abs(threshold), abs(threshold))
    } else if(selection == "lower") {
      threshold <- c(threshold, Inf)
    } else {
      threshold <- c(-Inf, threshold)
    }
  }

  lowerProb <- pnorm(threshold[1], sd = sqrt(contrastVar), lower.tail = TRUE)
  upperProb <- pnorm(threshold[2], sd = sqrt(contrastVar), lower.tail = FALSE)
  pthreshold <- lowerProb + upperProb

  # Validating selection --------------
  y <- as.numeric(y)
  testStat <- sum(y * contrast)
  if(testStat > threshold[1] & testStat < threshold[2]) {
    stop("Test statistic did not cross threshold!")
  }

  t2 <- confidence_level^2 * pthreshold

  # Computing MLE -------------------------
  if("mle" %in% estimate_type) {
    mleMu <- computeLinearMLE(y, sigma, contrast, threshold)
  } else {
    mleMu <- NULL
  }

  if(estimate_type[1] == "mle") {
    muhat <- mleMu
  }

  # Sampling from the global-null ---------------------
  if(is.null(nSamples)) {
    nSamples <- round(length(y) * 5 / confidence_level)
    nSamples <- max(nSamples, 100 * length(y))
  }

  if(any(c("global-null", "hybrid") %in% pvalue_type)) {
    nullMu <- rep(0, length(y))
    nullSample <- sampleLinearTest(nSamples, mu = rep(0, p), sigma, contrast, threshold)
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
    polyResult <- getPolyCI(y, sigma, contrast, threshold, confidence_level,
                            test = "linear")
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

  # Null Distribution Based CIs ---------------
  if(("global-null" %in% ci_type)) {
    nullCI <- linearRB(y, sigma, contrast, threshold,
                       confidence_level, rbIters = NULL,
                       variables = NULL, computeFull = TRUE)
  } else {
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
  if("switch" %in% ci_type) {
    if(verbose) print("Computing Switching Regime CIs!")
    switchCI <- getSwitchCI(y, sigma, contrast, threshold, pthreshold,
                            confidence_level, quadlam,
                            confidence_level^2, testStat,
                            hybridPval, trueHybrid, rbIters,
                            test = "linear")
  } else {
    switchCI <- NULL
  }

  if(ci_type[1] == "switch") {
    ci <- switchCI
  }

  # Computing hybrid CIs ---------------------
  if("hybrid" %in% ci_type) {
    if(verbose) print("Computing Hybrid CIs!")
    hybridCI <- getHybridCI(y, sigma, contrast, threshold, pthreshold, confidence_level,
                            hybridPval = hybridPval, trueHybrid, rbIters,
                            test = "linear")
  } else {
    hybridCI <- NULL
  }

  if(ci_type[1] == "hybrid") {
    ci <- hybridCI
  }

  # Wrapping up ------------------
  results <- list()
  results$muhat <- muhat
  results$ci <- ci
  results$pvalue <- pvalue

  results$naiveMu <- y
  results$mleMu <- mleMu
  results$nullSample <- nullSample
  results$nullPval <- nullPval
  results$polyPval <- polyPval
  results$polyCI <- polyCI
  results$nullCI <- nullCI
  results$naivePval <- naivePval
  results$naiveCI <- naiveCI
  results$switchCI <- switchCI
  results$hybridPval <- hybridPval
  results$hybridCI <- hybridCI

  results$contrast <- contrast
  results$sigma <- sigma
  results$pthreshold <- pthreshold
  results$threshold <- threshold
  results$confidence_level <- 1 - confidence_level
  results$switchTune <- t2
  results$estimate_type <- estimate_type
  results$pvalue_type <- pvalue_type
  results$ci_type <- ci_type
  results$testStat <- testStat
  results$trueHybrid <- trueHybrid
  results$rbIters <- rbIters

  class(results) <- "mvnLinear"

  return(results)
}

computeLinearMLE <- function(y, sigma, contrast, threshold) {
  naive <- sum(y * contrast)
  cSig <- as.numeric(t(contrast) %*% sigma)
  cSig <- cSig / as.numeric(cSig %*% contrast)
  interval <- c(-abs(naive), abs(naive))
  interval <- sort(interval)
  sd <- as.numeric(sqrt(t(contrast) %*% sigma %*% contrast))

  maximum <- optimize(f = computeLinearDens, interval = interval,
                      maximum = TRUE,
                      naive, y, cSig, sd, threshold, solve(sigma))$maximum
  return(y + (maximum - naive) * cSig)
}

computeLinearDens <- function(m, naive, y, cSig, sd, threshold, precision) {
  diff <- (m - naive) * cSig
  density <- -0.5 * t(diff) %*% precision %*% diff
  prob <- pnorm(threshold[1], sd = sd, mean = m) +
    pnorm(threshold[2], sd = sd, mean = m, lower.tail = FALSE)

  condDens <- density - prob
  return(condDens)
}





