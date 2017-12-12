#' Inference for Generalized Linear Models Selected via an Aggregate Test
#'
#' @description \code{psatGLM} is used to estimate regression models that were selected
#' based on a single aggregate test of the of form \deqn{\beta' K \beta > c > 0}
#' or \deqn{a'\beta < l or a'\beta > u}
#' where \eqn{\beta} is the MLE estimate for the regression coefficients.
#'
#' @param X the design matrix of the regression model without an
#' intercept.
#'
#' @param y the observed dependent variable.
#'
#' @param test For an aggregate test, the string `wald` or a semi-positive matrix
#' of dimension \code{ncol(X) X ncol{X}}.
#' For a linear test a vector of length \code{ncol(X)}.
#'
#' @param test_direction If a linear test is specified and no threshold is specififed
#' then this will be used to determine the threshold.
#'
#' @param family see \code{\link[stats]{glm}} for details.
#'
#' @param resid_sd method for estimating the covariance matrix. Re-fitted least squares
#' if `naive` or intercept model if `null`.
#'
#' @param threshold a threshold that the aggregate test must cross in order for the model
#' to be estimated. Should be a positive scalar for quadratic tests and a vector of size 2
#' for linear aggregate tests.
#'
#' @param pval_threshold if threshold is not specified, then a threshold will be set such that
#' the selection occures iff \code{aggregate_pval < pval_threshold}.
#'
#' @param estimate_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param pvalue_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param ci_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param confidence_level see \code{\link{mvnQuadratic}} for details.
#'
#' @param verbose whether the progress of the inference function should be reported.
#'
#' @param switch_tune see \code{\link{mvnQuadratic}} for details.
#'
#' @param nSamples see \code{\link{mvnQuadratic}} for details.
#'
#' @param verbose see \code{\link{mvnQuadratic}} for details.
#'
#' @param control an object of type \code{\link{psatControl}}.
#'
#' @details \code{psatGLM} is used to perform inference in generalized
#' linear models that were selected via a single aggregate test. This
#' function essentially works as a wrapper for \code{\link{mvnQuadratic}}
#' and \code{\link{mvnLinear}}. See the
#' documentation of \code{\link{mvnQuadratic}} for details regrading the
#' available inference methods.
#'
#' The function works by fitting a regression model using the \code{\link[stats]{glm}}
#' function. If \code{resid_sd} is set to "naive" then the covariance matrix is estimated
#' based on the initial model fit. If \code{resid_sd} is set to "null" then the covariance matrix
#' will be estimated based on an intercept model. Once a model has been estimated, an
#' aggregate test is performed based on the estimated covariance with the intercept removed
#' and the vector of regression coefficients is passed to \code{\link{mvnQuadratic}} or \code{\link{mvnLinear}}
#' without the intercept.
#'
#' @return an object of type \code{psatGLM}
#'
#' @seealso \code{\link{mvnQuadratic}}, \code{\link{getCI}},
#' \code{\link{getPval}}, \code{\link{plot.psatGLM}},
#' \code{\link{summary.psatGLM}}, \code{\link{coef.pastGLM}},
#' \code{\link{predict.psatGLM}}
psatGLM <- function(X, y, test = "wald",
                    test_direction = c("two-sided", "lower", "upper"),
                    family = "gaussian",
                    resid_sd = c("null", "naive"),
                    threshold = NULL, pval_threshold = 0.05,
                    estimate_type = c("mle", "naive"),
                    pvalue_type = c("hybrid", "polyhedral", "naive"),
                    ci_type = c("switch", "polyhedral", "naive"),
                    confidence_level = .95,
                    verbose = TRUE,
                    control = psatControl()){

  family <- family[1]
  if(any(family %in% c("quasi", "quasibinomial", "quasipoisson"))) {
    warning("Overdispersed families are not explicity supported!")
  }

  # What is the screening method? ------------
  if(is.matrix(test)) {
    if(all(dim(test) == ncol(X)))
    testType <- "quadratic"
  } else if(length(test) == ncol(X)) {
    testType <- "linear"
  } else if(is.character(test)) {
    if(test == "wald") {
      testType <- "quadratic"
    }
  } else {
    stop("Unable to determine test type. `test' variable must be one of: (1) A matrix of dimension ncol(X) x ncol(X) (for a general quadratic test). (2) The string `wald' (for Wald test). (3) A vector of length ncol(X).")
  }

  if(testType == "quadratic" & length(test_direction) == 1) {
    warning("`test_direction' is irrelevant for a quadratic test!")
  }

  if(testType == "quadratic" & length(threshold) > 1) {
    stop("For a quadratic test, `threshold' must be either a scalar or NULL!")
  }

  # Starting analysis ------------
  naivefit <- glm(y ~ X, family = family)
  naiveBeta <- coef(naivefit)[-1]

  resid_sd <- resid_sd[1]
  if(is.character(resid_sd)) {
    if(resid_sd == "null") {
      nullfit <- glm(y ~ 1, family = family)
      w <- vcov(nullfit) * length(y)
      sigma <- solve(t(X) %*% X) * as.numeric(w)
    } else if(resid_sd == "naive") {
      sigma <- vcov(naivefit)[-1, , drop = FALSE][, -1, drop = FALSE]
    } else {
      stop("Variance estimation method not supported!")
    }
  } else if(is.numeric(resid_sd) & resid_sd >0){
    sigma <- solve(t(X) %*% X) * resid_sd^2
  } else {
    stop("Invalid resid_sd value!")
  }

  if(testType == "quadratic") {
    mvnfit <- mvnQuadratic(naiveBeta, sigma, testMat = test,
                           threshold = threshold, pval_threshold = pval_threshold,
                           estimate_type = estimate_type,
                           pvalue_type = pvalue_type,
                           ci_type = ci_type,
                           confidence_level = confidence_level,
                           verbose = verbose, control = control)
  } else if(testType == "linear") {
    mvnfit <- mvnLinear(y = naiveBeta, sigma = sigma, contrast = test,
                        threshold = threshold, pval_threshold = pval_threshold,
                        selection = test_direction,
                        estimate_type = estimate_type,
                        pvalue_type = pvalue_type,
                        ci_type = ci_type,
                        confidence_level = confidence_level,
                        verbose = verbose, trueHybrid = TRUE)
  }

  results <- mvnfit
  results$y <- y
  results$X <- X
  results$naiveBeta <- coef(naivefit)

  if(!is.null(results$mleMu)) {
    beta <- results$mleMu
    offset <- as.numeric(X %*% beta)
    interceptfit <- glm(y ~ offset(offset), family = family)
    intercept <- coef(interceptfit)[1]
    interceptsd <- summary(interceptfit)$coefficients[1, 2]
    beta <- c(intercept, beta)
    results$mleBeta <- beta
  } else {
    interceptsd <- summary(naivefit)$coefficients[1, 2]
  }

  if("mle" %in% results$estimate_type) {
    results$betahat <- results$mleBeta
  } else {
    results$betahat <- results$naiveBeta
  }

  results$muhat <- NULL
  results$mleMu <- NULL
  results$naiveMu <- NULL
  results$family <- family
  results$interceptsd <- interceptsd
  class(results) <- c("psatGLM", class(results))
  return(results)
}
