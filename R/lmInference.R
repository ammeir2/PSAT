glmQuadratic <- function(X, y, testMat = "wald", family = "gaussian",
                         resid_sd = c("ysd", "naive"),
                         threshold = NULL, pval_threshold = 0.05,
                         estimate_type = c("mle", "naive"),
                         pvalue_type = c("hybrid", "polyhedral", "naive", "global-null"),
                         ci_type = c("switch", "polyhedral", "naive", "global-null"),
                         confidence_level = 0.05,
                         switchTune = c("sqrd", "half"),
                         nSamples = NULL, verbose = TRUE) {
  family <- family[1]
  if(any(family %in% c("quasi", "quasibinomial", "quasipoisson"))) {
    warning("Overdispersed families are not explicity supported!")
  }

  naivefit <- glm(y ~ X, family = family)
  naiveBeta <- coef(naivefit)[-1]
  sigma <- vcov(naivefit)[-1, , drop = FALSE][, -1, drop = FALSE]

  resid_sd <- resid_sd[1]
  if(is.character(resid_sd)) {
    if(resid_sd == "ysd") {
      ysd <- sd(y)
    } else if(resid_sd == "naive") {
      ysd <- sd(naivefit$residuals)
    } else {
      stop("Variance estimation method not supported!")
    }
  }

  if(family == "gaussian") {
    sigma <- sigma * ysd^2 / var(naivefit$residuals)
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
