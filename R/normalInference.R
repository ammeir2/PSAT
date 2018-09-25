#' Inference for Normal Means after Aggregate Testing
#'
#' @description \code{mvnQuadratic} is used to estimate a normal means model that was selected based
#' on a single quadratic aggregate test of the form:
#' \deqn{y' K y > c > 0}
#'
#' @param y the observed normal vector.
#'
#' @param sigma the covariance matrix of \code{y}.
#'
#' @param testMat the test matrix \eqn{K} used in the aggregate test
#'
#' @param threshold the threshold \eqn{c > 0} used in the aggregate test.
#'
#' @param pval_threshold the signficance level of the aggregate test.
#' Overrided by \code{threshold} if both are provided.
#' 
#' @param contrasts an optional matrix of contrasts to be tested: must have number of columns 
#' identical to the length of \code{y}. If left as \code{NULL}, the coorinates of \code{y} will be tested by default. 
#'
#' @param estimate_type the types of point estimates to compute and report. The first
#' estimator listed will be used as the default method.
#'
#' @param pvalue_type a vector of methods with which to compute the p-values. The first
#' method listed will be used as the default method.
#'
#' @param ci_type a vector of confidence interval computation methods to be used.
#'  The first method listed will be will be used as the default method.
#'
#' @param confidence_level the confidence level for constructing confidencei intervals.
#'
#' @param verbose whether to report on the progress of the computation.
#'
#' @param control an object of type \code{\link{psatControl}}.
#'
#' @details The function is used to perform inference for normal mean vectors
#' that were selected based on a single quadratic aggregate test. To be exact, suppose
#' that \eqn{y ~ N(\mu,\Sigma)} and that we are interested in estimating \eqn{\mu}
#' only if we can determine that \eqn{\mu\neq 0} using an aggregate test of the form:
#' \deqn{y' K y > c > 0} for some predetermined constant \eqn{c}. If \code{testMat} is set
#' to the default value of "wald", then \eqn{K = \Sigma^{-1}}. If wald test is used, it is
#' recommended to specify \code{testMat} as "wald" because this setting makes some of computations
#' more efficient. Otherwise, \code{testMat} must be a positive definite matrix of an
#' appropriate dimension.
#'
#' If \code{estimate_type} includes the string "mle" then \code{mvnQuadratic}
#' will compute the conditional maximum likelihood estimator for the mean vector,
#' which is typically a shrinkage estimator.  If \code{testMat = "wald"} then the
#' computation is performed via an efficient line-search method. Otherwise, 
#' the computation is performed via the Nelder-Mead method where the probability
#' of selection is approximated using the \code{\link[CompQuadForm]{liu}} function.
#' 
#' The \code{threshold} parameter specifies the constant \eqn{c>0} which is used
#' to threshold the aggregate test. It takes precedence over \code{pval_threshold} if both
#' are specified. We use the \code{\link[CompQuadForm]{liu}} function to compute the the
#' threshold if a non-Wald test is used.
#' 
#' \code{mvnQuadratic} offers several options for computing p-values. The "global-null"
#' method relies on comparing the magnitude of \eqn{y} to samples from the truncated
#' global-null distribution. This method is powerful when \eqn{\mu} is sparse and its
#' non-zero coordinates are not very large. The "polyhedral" method is exact when the
#' observed data is approximately normal and is quite robust to model misspecification.
#' It tends to be more powerful than the `global-null` method when the magnitude of
#' \eqn{\mu} is large. The "hybrid" method combines the strengths of the "global-null"
#' and "polyhedral" methods, possessing good power regardless of the sparsity or
#' magnitude of \eqn{\mu}. However it is less robust to the misspecification of the distribution
#' of \eqn{y} than the "polyhedral" method. The confidence interval methods are similar to the p-values ones,
#' with the Regime switching
#' confidence intervals ("switch") serving a simialr purpose as the "hybrid" method.
#'
#' @return An object of class \code{mvnQuadratic}.
#'
#' @seealso \code{\link{getCI}}, \code{\link{getPval}},
#' \code{\link{coef.mvnQuadratic}}, \code{\link{plot.mvnQuadratic}},
#' \code{\link{psatGLM}}.
mvnQuadratic <- function(y, sigma, testMat = "wald",
                         threshold = NULL, pval_threshold = 0.05,
                         contrasts = NULL,
                         estimate_type = c("mle", "naive"),
                         pvalue_type = c("hybrid", "polyhedral", "naive"),
                         ci_type = c("switch", "polyhedral", "naive"),
                         confidence_level = .95,
                         verbose = TRUE,
                         control = psatControl()) {
  # Getting control variables
  if(class(control) != "quadratic_control") {
    stop("control must be an object of class `quadratic_control'!")
  }

  switchTune <- control$switchTune
  nullMethod <- control$nullMethod
  nSamples <- control$nSamples
  sgdStep <- control$sgdStep
  nsteps <-  control$nsteps
  trueHybrid <- control$trueHybrid
  rbIters <- control$rbIters
  truncPmethod <- control$truncPmethod
  quadraticSampler <- control$quadraticSampler

  # Validating parameters --------------------
  if(!(nullMethod %in% c("zero-quantile", "RB"))) {
    stop("nullMethod must be either `zero-quantile' or `RB'!")
  }

  p <- length(y)
  if(!(any(estimate_type %in% c("mle", "naive")))) {
    stop("estimation method not supported!")
  }

  checkPvalues(pvalue_type)
  checkCI(ci_type)

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

  if(confidence_level <= 0 | confidence_level >= 1) {
    stop("confidence_level must be between 0 and 1!")
  }
  confidence_level <- 1 - confidence_level

  if(is.null(switchTune)) {
    switchTune <- confidence_level^2
  }
  
  # Contrasts ------
  if(is.null(contrasts)) {
    contrasts <- diag(length(y))
  }
  
  if(ncol(contrasts) != length(y) | 
     nrow(contrasts) < 1) {
    stop("Contrasts must be a matrix with length(y) columns! (and at least one row)")
  }

  # Setting threshold ----------------------------
  p <- nrow(sigma)
  pthreshold <- pval_threshold
  if(is.null(threshold)) {
    if(pthreshold < 0 | pthreshold > 1) {
      stop("If a threshold is not provided, then pval_threshold must be between 0 and 1!")
    } else {
      threshold <- getQuadraticThreshold(pthreshold, quadlam)
    }
  }

  if(is.null(pval_threshold)) {
    pthreshold <- CompQuadForm::liu(threshold, quadlam)
  }

  y <- as.numeric(y)
  testStat <- as.numeric(t(y) %*% testMat %*% y)
  if(testStat < threshold) {
    warning("Test statistic is below the computed threshold, verify that the model was selected!")
    warning("Assuming that the effective threshold is the observed value.")
    threshold <- testStat - testStat/10^-4
    pthreshold <- CompQuadForm::liu(threshold, quadlam)
  }

  t2 <- switchTune * pthreshold

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
      if(control$optimMethod == "SGD") {
        if(any(quadlam < 10e-8)) {
          control$optimMethod <- "Nelder-Mead"
          warning("Test matrix nearly singualr: using Nelder-Mead for MLE computation")
        } else {
          sgdfit <- quadraticSGD(y, sigma, invcov, testMat, threshold,
                                 stepRate = 0.75, stepCoef = sgdStep,
                                 delay = 20, sgdSteps = nsteps, assumeCovergence = 800,
                                 mhIters = 40)
          mleMu <- sgdfit$mle
          solutionPath <- sgdfit$solutionPath
        }
      } 
      
      if(control$optimMethod == "Nelder-Mead") {
        mleMu <- quadraticNM(y, sigma, invcov, testMat, threshold)
        solutionPath <- NULL
      }
      
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
    nSamples <- round(length(y) * 5 / confidence_level)
    nSamples <- max(nSamples, 30 * length(y))
  }

  if(any(c("global-null", "hybrid") %in% pvalue_type) |
     (nullMethod == "zero-quantile")) {
    if(verbose) print(paste("Sampling from null distribution! (", nSamples, " samples)", sep = ""))
    nullMu <- rep(0, length(y))
    
    if(quadraticSampler == "PSAT") {
      if(any(quadlam < 10e-8)) {
        warning("Test matrix nearly singular: using tmg sampler instead of PSAT sampler!")
        quadraticSampler <- "tmg"
      } else {
        nullSample <- sampleQuadraticConstraint(nullMu, as.matrix(sigma),
                                                init = y,
                                                threshold, testMat,
                                                sampSize = nSamples,
                                                burnin = 1000,
                                                trim = 50)
      }
    }
    
    if(quadraticSampler == "tmg") {
      nullSample <- rtmg(n = nSamples * 10,
                         M = invcov,
                         r = as.numeric(invcov %*% nullMu),
                         initial = y,
                         q = list(a = list(A = testMat, B = rep(0, p), C = -threshold)),
                         burn.in = 200)
    }
    nullPval <- numeric(nrow(contrasts))
    
    for(i in 1:nrow(contrasts)) {
      contt <- sum(contrasts[i, ] * y)
      contsamp <- as.numeric(nullSample %*% contrasts[i, ])
      nullPval[i] <- 2 * min(mean(contt < contsamp), mean(contt > contsamp))
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
    polyResult <- suppressWarnings(getPolyCI(y, sigma, testMat, threshold, confidence_level,
                                             contrasts = contrasts,
                                             truncPmethod = truncPmethod))
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

  # MLE based inference (removed) -----------
  mleSample <- NULL
  mlePval <- NULL
  mleCI <- NULL
  mleDist <- NULL

  # Null Distribution Based CIs ---------------
  if(("global-null" %in% ci_type)) {
    if(nullMethod == "zero-quantile") {
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
      nullCI <- quadraticRB(y, sigma, testMat, threshold,
                            variables = contrasts,
                            confidence_level, computeFull = TRUE,
                            rbIters = rbIters)
      nullDist <- NULL
    }
  } else {
    nullDist <- NULL
    nullCI <- NULL
  }

  if(ci_type[1] == "global_null") {
    ci <- nullCI
  }

  # Naive pvalues and CIs -----------------------
  naiveCI <- getNaiveCI(y, sigma, confidence_level, contrasts, TRUE)
  naivePval <- naiveCI$pval
  naiveCI <- naiveCI$ci
  if(ci_type[1] == "naive") {
    ci <- naiveCI
  }
  if(pvalue_type[1] == "naive") {
    pvalue <- naivePval
  }
  
  # Regime switching CIs -------------------------
  if("switch" %in% ci_type) {
    if(verbose) print("Computing Switching Regime CIs!")
    switchCI <- getSwitchCI(y, sigma, testMat, threshold, pthreshold,
                            confidence_level, contrasts = contrasts,
                            quadlam,
                            t2 = pthreshold * switchTune, testStat,
                            hybridPval, trueHybrid, rbIters)
  } else {
    switchCI <- NULL
  }

  if(ci_type[1] == "switch") {
    ci <- switchCI
  }

  # Computing hybrid CIs ---------------------
  if("hybrid" %in% ci_type) {
    if(verbose) print("Computing Hybrid CIs!")
    hybridCI <- getHybridCI(y, sigma, testMat, threshold, pthreshold, confidence_level,
                            hybridPval = hybridPval, trueHybrid, rbIters)
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
  results$contrasts <- contrasts

  results$y <- y
  results$naiveContrast <- as.numeric(contrasts %*% y)
  results$mleMu <- mleMu
  if(!is.null(mleMu)) {
    results$mleContrast <- as.numeric(contrasts %*% mleMu)
  } 
  results$solutionPath <- solutionPath
  results$nullSample <- nullSample
  results$nullDist <- nullDist
  results$nullPval <- pmax(nullPval, naivePval)
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
  results$hybridCI <- hybridCI

  results$testMat <- testMat
  results$sigma <- sigma
  results$invcov <- invcov
  results$pthreshold <- pthreshold
  results$threshold <- threshold
  results$confidence_level <- 1 - confidence_level
  results$switchTune <- switchTune
  results$estimate_type <- estimate_type
  results$pvalue_type <- pvalue_type
  results$ci_type <- ci_type
  results$nullMethod <- nullMethod
  results$testStat <- testStat
  results$quadlam <- quadlam
  results$trueHybrid <- trueHybrid
  results$rbIters <- rbIters
  results$testType <- "quadratic"

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

#' Creates a list of parameters for use with PSAT inference functions
#' 
#' Creates a list with additional parameters for use with 
#' \code{\link{mvnQuadratic}}, \code{\link[psat]{mvnLinear}}, and \code{\link[PSAT]{psatGLM}}.
#' 
#' @param switchTune tuning parameter for regime switching confidence intervals, 
#' should be a number between 0 and confidence_level.
#' 
#' @param nullMethod method for compute global-null confidence intervals. Robins-Monroe (RB) should be used.
#' 
#' @param nSamples number of samples to be taken from the null distribution.
#' 
#' @param sgdStep number of stochastic gradient steps to take, only applicable for 
#' non-wald quadratic tests when \code{optimMethod} "SGD" is used.
#' 
#' @param trueHybrid whether Robins-Monroe computation should be performed for computing 
#' confidence intervals for all confidence intervals or only when they may improve power.
#' 
#' @param rbIters number of steps to take when computing confidence intervals with the Robins-Monroe
#' procedure.
#' 
#' @param optimMethod optimization method to be used when computing the conditional MLE in
#' in inference after testing with a non-wald aggregate test. Nelder-Mead as implemented in the 
#' \code{\link[stats]{optim}} function.
#' 
#' @param truncPmethod the type of test to use when computing the polyhedral p-values,
#' options are either the UMPU test or a symmetric test. 
#' 
#' @param quadraticSampler which quadratic sampler to use? Choices are either the 
#' Hamiltionian Montel-Carlo method implemented in \code{\link[tmg]{rtmg}},
#' or a Gibbs sampler implemented in this package. 
psatControl <- function(switchTune = NULL,
                        nullMethod = c("RB", "zero-quantile"),
                        nSamples = NULL,
                        sgdStep = NULL,
                        nsteps = 1000,
                        trueHybrid = FALSE,
                        rbIters = NULL,
                        optimMethod = c("Nelder-Mead", "SGD"),
                        truncPmethod = c("UMPU", "symmetric"),
                        quadraticSampler = c("tmg", "PSAT")) {
  control <- list()
  control$switchTune <- switchTune
  control$nullMethod <- nullMethod[1]
  control$nSamples <- nSamples
  control$sgdStep <- sgdStep
  control$nsteps <- nsteps
  control$trueHybrid <- trueHybrid
  control$rbIters <- rbIters
  control$optimMethod <- optimMethod[1]
  control$truncPmethod <- truncPmethod[1]
  control$quadraticSampler <- quadraticSampler[1]
  class(control) <- "quadratic_control"
  return(control)
}

checkPvalues <- function(pvalvec) {
  if(!(any(pvalvec %in% c("polyhedral", "hybrid", "global-null", "mle", "naive")))) {
    stop("pvalue_type not supported!")
  }
}

checkCI <- function(civec) {
  if(!(any(civec %in% c("polyhedral", "switch", "naive", "mle", "hybrid")))) {
    stop("ci_type not supported!")
  }
}


