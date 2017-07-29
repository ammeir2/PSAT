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
    llim <- theta
    tempPval <- ptruncNorm(llim, theta, sqrt(etaSigma), lower, upper)
    while(tempPval < 1 - alpha / 2) {
      llim <- llim - sqrt(etaSigma) * 1
      tempPval <- ptruncNorm(llim, theta, sqrt(etaSigma), lower, upper)
    }

    ulim <- theta
    tempPval <- ptruncNorm(ulim, theta, sqrt(etaSigma), lower, upper)
    while(tempPval > alpha / 2) {
      ulim <- ulim + sqrt(etaSigma) * 0.1
      tempPval <- ptruncNorm(ulim, theta, sqrt(etaSigma), lower, upper)
    }

    uci <- uniroot(f = function(x) ptruncNorm(x, theta, sqrt(etaSigma), lower, upper) - alpha / 2,
                   interval = c(llim, ulim))$root
    lci <- uniroot(f = function(x) ptruncNorm(x, theta, sqrt(etaSigma), lower, upper) - (1 - alpha / 2),
                   interval = c(llim, ulim))$root
    ci <- c(lci, uci)
  } else {
    ci <- NULL
  }

  return(list(pval = etaPval, ci = ci))
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
