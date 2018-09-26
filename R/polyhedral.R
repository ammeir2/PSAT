polyhedral.workhorse <- function(eta, u, sigma, testMat, threshold,
                                 computeCI = FALSE, alpha = 0.05,
                                 truncPmethod = "symmetric") {
  # Computing truncation ----------------
  c <- as.numeric(as.numeric(t(eta) %*% sigma %*% eta)^-1 * sigma %*% eta)
  v <- as.numeric(u - c %*% t(eta) %*% u)
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
  if(truncPmethod == "symmetric" | delta < 0) {
    etaPval <- ptruncNorm(0, theta, sqrt(etaSigma), lower, upper)
    etaPval <- 2 * min(etaPval, 1 - etaPval)
  } else {
    etaPval <- computeUMPU(theta, lower, upper, 0, etaSigma)
  }

  # Computing CI ----------------
  if(computeCI) {#& truncPmethod == "symmetric") {
    ci <- suppressWarnings(findPolyCIlimits(theta, etaSigma, lower, upper, alpha))
  } else if(computeCI & truncPmethod == "UMPU" & FALSE) {
    ci <- UMAUsearch(theta, lower, upper, etaSigma, alpha)
  } else {
    ci <- NULL
  }

  return(list(pval = as.numeric(etaPval), ci = ci, noCorrection = delta < 0))
}

linearPolyhedral <- function(eta, u, sigma, a, threshold,
                             computeCI = FALSE, alpha = 0.05,
                             truncPmethod = "symmetric") {
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
  if(truncPmethod == "symmetric") {
    etaPval <- ptruncNorm(0, theta, sqrt(etaSigma), lower, upper)
    etaPval <- 2 * min(etaPval, 1 - etaPval)
  } else {
    etaPval <- computeUMPU(theta, lower, upper, 0, etaSigma)
  }

  # Computing CI ----------------
  if(computeCI) {
    ci <- suppressWarnings(findPolyCIlimits(theta, etaSigma, lower, upper, alpha))
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

# A function for computing 
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

dtruncNorm <- function(x, mu, sd, l, u)  {
  if(x > l & x < u) {
    return(0)
  } else {
    return(dnorm(x, mu, sd) / (pnorm(l, mu, sd) + 1 - pnorm(u, mu, sd)))
  }
}

# Function for computing the conditional expectation for a test-statistic
computeCondExp <- function(lower, upper, mu = 0, sd = 1, side = "upper") {
  ppos <- pnorm(upper, mean = mu, sd = sd, lower.tail = FALSE, log.p = TRUE) 
  pneg <- pnorm(lower, mean = mu, sd = sd, lower.tail = TRUE, log.p = TRUE)
  ppos <- 1 / (1 + exp(pneg - ppos))
  pneg <- 1 - ppos
  if(side == "upper") {
    return(etruncnorm(a = upper, b = Inf, mean = mu, sd = sd) * ppos)
  } else if(side == "lower") {
    return(etruncnorm(a = -Inf, b = lower, mean = 0, sd = sd) * pneg)
  } else {
    l <- etruncnorm(a = -Inf, b = lower, mean = 0, sd = sd) * pneg
    u <- etruncnorm(a = upper, b = Inf, mean = mu, sd = sd) * ppos
    return(l + u)
  }
}

# c1 and c2 are test limits
# lower and upper are polyherdal truncations (through conditioning on selection event and W)
computeTestExp <- function(c1, c2, lower, upper, mu, sd) {
  # c1 <- c2 - exp(c1)
  if(c1 < upper & c1 > lower) {
    c1 <- lower
  }
  
  if(c2 < lower) {
    pupper <- pnorm(upper, mean = mu, sd = sd, lower.tail = FALSE)
    pmid <- pnorm(lower, mean = mu, sd = sd) - pnorm(c2, mean = mu, sd = sd)
    midExp <- etruncnorm(c2, lower, mean = mu, sd = sd)
    upperExp <- (etruncnorm(upper, Inf, mean = mu, sd = sd) * pupper + midExp * pmid) / (pupper + pmid)
    pupper <- pupper + pmid
  } else {
    pupper <- pnorm(c2, mean = mu, sd = sd, lower.tail = FALSE)
    upperExp <- etruncnorm(c2, Inf, mean = mu, sd = sd)
  }
  
  if(c1 > upper) {
    plower <- pnorm(lower, mean = mu, sd = sd, lower.tail = FALSE)
    pmid <- pnorm(c1, mean = mu, sd = sd) - pnorm(upper, mean = mu, sd = sd)
    midExp <- etruncnorm(upper, c1, mean = mu, sd = sd)
    lowerExp <- (etruncnorm(-Inf, lower, mean = mu, sd = sd) * plower + midExp * pmid) / (plower + pmid)
    plower <- plower + pmid
  } else {
    plower <- pnorm(c1, mean = mu, sd = sd, lower.tail = TRUE)
    lowerExp <- etruncnorm(-Inf, b = c1, mean = mu, sd = sd)
  }
  
  expectation <- (plower * lowerExp + pupper * upperExp) / (plower + pupper)
  pselect <- pnorm(lower, mean = mu, sd = sd) + pnorm(upper, mean = mu, sd = sd, lower.tail = FALSE)
  return(list(exp = expectation, level = (plower + pupper) / pselect))
}

# Optima Test equation
optimTestEquation <- function(c1, c2, tlower, tupper, mu, sd, truncExp) {
  levelExp <- computeTestExp(c1, c2, tlower, tupper, mu, sd) 
  # return((levelExp$exp - levelExp$level * truncExp)^2)
  return(levelExp$exp - levelExp$level * truncExp)
}

# A Function for computing Fithian's optimal test
computeUMPU <- function(theta, lower, upper, mu, etaSigma) {
  if(theta <= mu) {
    c2 <- -theta
    tlower <- -upper
    tupper <- -lower
  } else {
    tlower <- lower
    tupper <- upper
    c2 <- theta
  }
  
  # c1 <- nlm(f = optimTestEquation, p = log(2 * (c2 - mu)),
  #           c2 = c2, tlower = tlower, tupper = tupper,
  #           mu = mu, sd = sqrt(etaSigma), truncExp = condExp)$estimate
  # testLevel <- computeTestExp(c1, c2, tlower, tupper, mu, sqrt(etaSigma))
  # c1 <- c2 - exp(c1)
  
  condExp <- computeCondExp(tlower, tupper, mu, sqrt(etaSigma), side = "both")
  c1 <- UMPUlowerValSearch(c2, tlower, tupper, mu, sqrt(etaSigma), condExp)
  testLevel <- computeTestExp(c1, c2, tlower, tupper, mu, sqrt(etaSigma))
  return(testLevel$level)
}

UMPUlowerValSearch <- function(c2, l, u, mu, sd, truncExp) {
  c1 <- mu
  fval <- optimTestEquation(mu, c2, l, u, mu, sd, truncExp)
  # if(fval < 0) {
  #   return(mu - (c2 - mu))
  # }
  nval <- fval
  while(sign(fval) == sign(nval)) {
    c1 <- c1 - sd
    nval <- optimTestEquation(c1, c2, l, u, mu, sd, truncExp)
    if(c1 < mu - 10 * sd) {
      c1 <- min(l, mu - (c2 - mu))
      return(c1)
    }
  }
  
  uniresult <- uniroot(f = optimTestEquation, interval = c(c1, mu),
                       c2 = c2, tlower = l, tupper = u, mu = mu, sd = sd,
                       truncExp = truncExp)
  return(uniresult$root)
}

UMAUsearch <- function(theta, lower, upper, etaSigma, alpha) {
  maxSteps <- 5
  ulimit <- theta + 0.01 * sqrt(etaSigma)
  upval <- computeUMPU(theta, lower, upper, ulimit, etaSigma)
  ucount <- 0
  while(upval < 1 - alpha / 2) {
    ulimit <- ulimit + sqrt(etaSigma)
    upval <- computeUMPU(theta, lower, upper, ulimit, etaSigma)
    ucount <- ucount + 1
    if(ucount > maxSteps) {
      break
    }
  }
  
  llimit <- theta - 0.01 * sqrt(etaSigma) 
  lpval <- computeUMPU(theta, lower, upper, llimit, etaSigma)
  lcount <- 0
  while(lpval > alpha / 2) {
    llimit <- llimit - sqrt(etaSigma)
    lpval <- computeUMPU(theta, lower, upper, llimit, etaSigma)
    lcount <- lcount + 1
    if(lcount > maxSteps) {
      break
    }
  }
  
  if(ucount <= maxSteps) {
    uci <- uniroot(function(m) computeUMPU(theta, lower, upper, m, etaSigma) - (1 - alpha / 2),
                   interval = c(llimit, ulimit))$root
  } else {
    uci <- ulimit
  }
  
  if(lcount <= maxSteps) {
    lci <- uniroot(function(m) computeUMPU(theta, lower, upper, m, etaSigma) - alpha / 2,
                   interval = c(llimit, ulimit))$root
  } else {
    lci <- llimit
  }
  
  return(c(lci, uci))
}


