#' Inference for Normal Means after Aggregate Testing with a Linear Test
#'
#' @description \code{mvnLinear} is used to estimate a normal means model that was selected based
#' on a single linear aggregate test of the form:
#' \deqn{a'y > u or a'y < l,} 
#' \deqn{l < u.}
#'
#' @param y the observed normal vector.
#'
#' @param sigma the covariance matrix of \code{y}.
#'
#' @param contrast the test contrast \eqn{a} of size \code{length(y)} used in the aggregate test
#'
#' @param threshold the threshold used in the aggregate test, either a vector of size two for the lower and upper
#' thresholds \eqn{u, l}, or a single number.
#'
#' @param pval_threshold the signficance level of the aggregate test.
#' Overrided by \code{threshold} if both are provided.
#' 
#' @param test_direction whether the linear test is one-sided or two-sided. Will be used if the provided threshold 
#' is a scalar, if lower then the tests will be \eqn{a'y < threshold} and if upper then the test will be \eqn{a'y > threshold}.
#'
#' @param estimate_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param pvalue_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param ci_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param confidence_level the confidence level for constructing confidencei intervals.
#'
#' @param verbose whether to report on the progress of the computation.
#' 
#' @param control an object of type \code{\link{psatControl}}.
#'
#' @details The function is used to perform inference for normal mean vectors
#' that were selected based on a single linear aggregate test. To be exact, suppose
#' that \eqn{y ~ N(\mu,\Sigma)} and that we are interested in estimating \eqn{\mu}
#' only if we can determine that \eqn{\mu != 0} using an aggregate test of the form:
#' \eqn{a'y <l} or \eqn{a'y > u} for some predetermined constants \eqn{a, l, u}.
#'
#' The \code{threshold} parameter specifies the constants \eqn{l<u} which are used
#' to threshold the aggregate test. If only a single number is provided, then the threshold
#' will be set according to test_direction: 
#' \itemize{
#' \item lower: a'y < threshold
#' \item upper: a'y > threshold
#' \item two-sided a'y < -threshold, or a'y > threshold
#' }
#' The \code{threshold} parameter takes precedence over \code{pval_threshold} if both
#' are specified. 
#'
#' See \code{\link{mvnQuadratic}} for details regarding the available inference methods.
#'
#' @return An object of class \code{mvnLinear}.
#'
#' @seealso \code{\link{mvnQuadratic}}, \code{\link{psatGLM}}, \code{\link{getCI}},
#' \code{\link{getPval}}.
mvnLinear <- function(y, sigma, contrast,
                      threshold = NULL, pval_threshold = 0.05,
                      test_direction = c("two-sided", "lower", "upper"),
                      estimate_type = c("mle", "naive"),
                      pvalue_type = c("hybrid", "polyhedral", "naive"),
                      ci_type = c("switch", "polyhedral", "naive"),
                      confidence_level = .95,
                      verbose = TRUE,
                      control = psatControl()) {

  # Getting control parameters ------
  switchTune <- control$switchTune
  nullMethod <- control$nullMethod
  nSamples <- control$nSamples
  trueHybrid <- control$trueHybrid
  rbIters <- control$rbIters

  # Checking input ----
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
  if(!any(test_direction %in% c("two-sided", "lower", "upper"))) {
     if(missing(threshold)) {
       stop("test_direction must be one of `two-sided', `lower' or `upper'!
         (or a threshold of length 2 must be provided)")
     } else if(length(threshold) < 2) {
       stop("test_direction must be one of `two-sided', `lower' or `upper'!
         (or a threshold of length 2 must be provided)")
     }
  }

  p <- length(y)
  contrastVar <- as.numeric(t(contrast) %*% sigma %*% contrast)
  test_direction <- test_direction[1]
  if(missing(threshold) | is.null(threshold)) {
    if(pthreshold < 0 | pthreshold > 1) {
      stop("If a threshold is not provided, then pval_threshold must be between 0 and 1!")
    } else {
      if(test_direction == "two-sided") {
        threshold <- qnorm(1 - pval_threshold/2, sd = sqrt(contrastVar))
        threshold <- c(-threshold, threshold)
      } else if(test_direction == "lower") {
        threshold <- qnorm(pval_threshold, sd = sqrt(contrastVar))
        threshold <- c(threshold, Inf)
      } else if(test_direction == "upper") {
        threshold <- qnorm(1 - pval_threshold, sd = sqrt(contrastVar))
        threshold <- c(-Inf, threshold)
      } else {
        stop("If threshold is not provided then test_direction must be one of `two-sided'
             `lower' or `upper'!")
      }
    }
  }

  if(length(threshold) > 2) {
    stop("The length of threshold must be either 1 or 2!")
  } else if(length(threshold) == 1) {
    if(test_direction == "two-sided") {
      warning(paste("Threshold of length one provided for a two sided test_direction rule.
                    test_direction rule assumed to be contrast < -abs(threshold) or contrast > abs(threshold)"))
      threshold <- c(-abs(threshold), abs(threshold))
    } else if(test_direction == "lower") {
      threshold <- c(threshold, Inf)
    } else {
      threshold <- c(-Inf, threshold)
    }
  }

  lowerProb <- pnorm(threshold[1], sd = sqrt(contrastVar), lower.tail = TRUE)
  upperProb <- pnorm(threshold[2], sd = sqrt(contrastVar), lower.tail = FALSE)
  pthreshold <- lowerProb + upperProb

  # Validating test_direction --------------
  y <- as.numeric(y)
  testStat <- sum(y * contrast)
  if(testStat > threshold[1] & testStat < threshold[2]) {
    stop("Test statistic did not cross threshold!")
  }

  # Regime switching CI tuning.
  if(is.null(switchTune)) {
    t2 <- confidence_level^2 * pthreshold
  } else {
    t2 <- switchTune * pthreshold
  }

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
                       confidence_level, rbIters = rbIters,
                       variables = NULL, computeFull = trueHybrid)
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
                            hybridPval = hybridPval, computeFull = trueHybrid, 
                            rbIters = rbIters,
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
  results$nullPval <- pmax(nullPval, naivePval)
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
  results$testType <- "linear"

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





