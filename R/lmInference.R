#' Inference for Generalized Linear Models Selected via an Aggregate Test
#'
#' @description \code{glmQuadratic} is used to estimate regression models that were selected
#' based on a single aggregate test of the of form \deqn{\beta' K \beta > c > 0}
#' where
#' \eqn{\hat\beta} is the MLE estimate for the regression coefficients.
#'
#' @param X the design matrix of the regression model without an
#' intercept.
#'
#' @param y the observed dependent variable.
#'
#' @param testMat see \code{\link{mvnQuadratic}} for details.
#'
#' @param family see \code{\link[stats]{glm}} for details.
#'
#' @param resid_sd method for estimating the covariance matrix.
#'
#' @param threshold see \code{\link{mvnQuadratic}} for details.
#'
#' @param pval_threshold see \code{\link{mvnQuadratic}} for details.
#'
#' @param estimate_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param pvalue_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param ci_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param confidence_level see \code{\link{mvnQuadratic}} for details.
#'
#' @param switch_tune see \code{\link{mvnQuadratic}} for details.
#'
#' @param nSamples see \code{\link{mvnQuadratic}} for details.
#'
#' @param verbose see \code{\link{mvnQuadratic}} for details.
#'
#' @details \code{glmQuadratic} is used to perform inference in generalized
#' linear models that were selected via a single aggregate test. This
#' function essentially works as a wrapper for \code{\link{mvnQuadratic}}. See the
#' documentaiton of \code{\link{mvnQuadratic}} for details regrading the inference
#' methods.
#'
#' The function works by fitting a regression model using the \code{\link[stats]{glm}}
#' function. The covariance matrix is estimated based on the estimated model if
#' \code{resid_sd} is set to "naive" or based on an intercept only model if
#' \code{resid_sd} is set to "null". Once a model has been estimated, a quadratic test
#' is performed based on the estimated covariance with the intercept removed
#' and the vector of regression coefficients is passed to \code{\link{mvnQuadratic}}
#' without the intercept.
#'
#' @return an object of type \code{glmQuadrtic}
#'
#' @seealso \code{\link{mvnQuadratic}}, \code{\link{getCI}},
#' \code{\link{getPval}}, \code{\link{plot.glmQuadratic}},
#' \code{\link{summary.glmQuadratic}}, \code{\link{coef.glmQuadratic}},
#' \code{\link{predict.glmQuadratic}}
glmQuadratic <- function(X, y, testMat = "wald", family = "gaussian",
                         resid_sd = c("null", "naive"),
                         threshold = NULL, pval_threshold = 0.05,
                         estimate_type = c("mle", "naive"),
                         pvalue_type = c("hybrid", "polyhedral", "naive", "global-null"),
                         ci_type = c("switch", "polyhedral", "naive", "global-null"),
                         confidence_level = .95,
                         switchTune = c("sqrd", "half"),
                         nSamples = NULL, verbose = TRUE) {
  family <- family[1]
  if(any(family %in% c("quasi", "quasibinomial", "quasipoisson"))) {
    warning("Overdispersed families are not explicity supported!")
  }

  naivefit <- glm(y ~ X, family = family)
  naiveBeta <- coef(naivefit)[-1]

  resid_sd <- resid_sd[1]
  if(is.character(resid_sd)) {
    if(resid_sd == "null") {
      nullfit <- glm(y ~ 1, family = family)
      w <- vcov(nullfit) * length(y)
      sigma <- solve(t(X) %*% X) * w
    } else if(resid_sd == "naive") {
      sigma <- vcov(naivefit)[-1, , drop = FALSE][, -1, drop = FALSE]
    } else {
      stop("Variance estimation method not supported!")
    }
  }


  mvnfit <- mvnQuadratic(naiveBeta, sigma, testMat = testMat,
                         threshold = threshold, pval_threshold = pval_threshold,
                         estimate_type = estimate_type,
                         pvalue_type = pvalue_type,
                         ci_type = ci_type,
                         confidence_level = confidence_level,
                         switchTune = switchTune,
                         nSamples = nSamples, verbose = verbose)
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
  class(results) <- "glmQuadratic"
  return(results)
}
