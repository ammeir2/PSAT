polyhedral.workhorse <- function(eta, u, sigma, testMat, threshold,
                                 computeCI = FALSE, alpha = 0.05) {
  # Computing truncation ----------------
  c <- as.numeric(t(eta) %*% sigma %*% eta)^-1 * sigma %*% eta
  v <- u - c %*% t(eta) %*% u
  vKc <- as.numeric(t(v) %*% testMat %*% c)
  cKc <- as.numeric(t(c) %*% testMat %*% c)
  vKv <- as.numeric(t(v) %*% testMat %*% v)

  delta <- 4 * (vKc)^2 - 4 * cKc * (vKv - threshold)
  if(delta >= 0) {
    upper <- (-2 * vKc + sqrt(delta)) / (2 * cKc)
    lower <- (-2 * vKc - sqrt(delta)) / (2 * cKc)
  } else {
    upper <- 0
    lower <- 0
  }

  # Computing p-value --------------
  theta <- as.numeric(t(eta) %*% u)
  etaSigma <- as.numeric(t(eta) %*% sigma %*% eta)
  etaPval <- ptruncNorm(0, theta, sqrt(etaSigma), lower, upper)
  etaPval <- 2 * min(etaPval, 1 - etaPval)

  # Computing CI ----------------
  if(computeCI) {
    ci <- findPolyCIlimits(theta, etaSigma, lower, upper, alpha)
  } else {
    ci <- NULL
  }

  return(list(pval = etaPval, ci = ci, noCorrection = delta < 0))
}

linearPolyhedral <- function(eta, u, sigma, a, threshold,
                             computeCI = FALSE, alpha = 0.05) {
  p <- length(u)
  cc <- as.numeric(t(eta) %*% sigma %*% eta)^-1 * sigma %*% eta
  theta <- sum(eta * u)
  w <- u - (cc * theta)
  lower <- as.numeric((threshold[1] - t(a) %*% w) / (t(a) %*% cc))
  upper <- as.numeric((threshold[2] - t(a) %*% w) / (t(a) %*% cc))
  if(lower > upper) {
    tempUpper <- upper
    upper <- lower
    lower <- tempUpper
  }

  etaSigma <- as.numeric(t(eta) %*% sigma %*% eta)
  etaPval <- ptruncNorm(0, theta, sqrt(etaSigma), lower, upper)
  etaPval <- 2 * min(etaPval, 1 - etaPval)

  # Computing CI ----------------
  if(computeCI) {
    ci <- findPolyCIlimits(theta, etaSigma, lower, upper, alpha)
  } else {
    ci <- NULL
  }

  return(list(pval = etaPval, ci = ci))
}

findPolyCIlimits <- function(theta, etaSigma, lower, upper, alpha) {
  llim <- theta
  ltempPval <- ptruncNorm(llim, theta, sqrt(etaSigma), lower, upper)
  while(ltempPval < 1 - alpha / 2) {
    llim <- llim - sqrt(etaSigma) * 0.1
    ltempPval <- ptruncNorm(llim, theta, sqrt(etaSigma), lower, upper)
    if(is.nan(ltempPval)) {
      llim <- llim + sqrt(etaSigma) * 0.1
      break
    }
  }

  ulim <- theta
  utempPval <- ptruncNorm(ulim, theta, sqrt(etaSigma), lower, upper)
  while(utempPval > alpha / 2) {
    ulim <- ulim + sqrt(etaSigma) * 0.1
    utempPval <- ptruncNorm(ulim, theta, sqrt(etaSigma), lower, upper)
    if(is.nan(utempPval)) {
      ulim <- ulim - sqrt(etaSigma) * 0.1
      break
    }
  }

  if(is.nan(utempPval)) {
    uci <- ulim
  } else {
    uci <- NULL
    capture.output(invisible(try(uci <- uniroot(f = function(x) ptruncNorm(x, theta, sqrt(etaSigma), lower, upper) - alpha / 2,
                   interval = c(llim, ulim))$root)))
    if(is.null(uci)){
      uci <- ulim
    }
  }

  if(is.nan(ltempPval)) {
    lci <- llim
  } else {
    lci <- NULL
    capture.output(invisible(try(lci <- uniroot(f = function(x) ptruncNorm(x, theta, sqrt(etaSigma), lower, upper) - (1 - alpha / 2),
                   interval = c(llim, ulim))$root)))
    if(is.null(lci)) {
      lci <- llim
    }
  }

  ci <- c(lci, uci)
  return(ci)
}

ptruncNorm <- function(mu, x, sd, l, u) {
  if(l == u) {
    return(pnorm(x, mu, sd))
  }
  prob <- 1 - pnorm(u, mu, sd) + pnorm(l, mu, sd)
  if(x < l) {
    return(pnorm(x, mu, sd) / prob)
  } else {
    return((pnorm(l, mu, sd) + pnorm(x, mu, sd) - pnorm(u, mu, sd)) / prob)
  }
}
