sqrtMat <- function(m) {
  eig <- eigen(m)
  sqrtm <- eig$vectors %*% diag(sqrt(abs(eig$values))) %*% t(eig$vectors)
  return(sqrtm)
}

sampleLinearTest <- function(nsamp, mu, sigma, contrast, threshold) {
  p <- length(mu)
  m <- sum(mu * contrast)
  sd <- as.numeric(sqrt(t(contrast) %*% sigma %*% contrast))
  negProb <- pnorm(threshold[1], mean = m, sd = sd, log.p = TRUE)
  posProb <- pnorm(threshold[2], mean = m, sd = sd, lower.tail = FALSE, log.p = TRUE)
  negProb <- 1 / (1 + exp(posProb - negProb))
  signs <- 1 - 2 * rbinom(nsamp, 1, negProb)
  contSamp <- numeric(nsamp)
  lower <- rep(-Inf, nsamp)
  lower[signs == 1] <- threshold[2]
  upper <- rep(threshold[1], nsamp)
  upper[signs == 1] <- Inf
  contSamp <- rtruncnorm(nsamp, a = lower, b = upper,
                         mean = m, sd = sd)

  cSig <- sigma %*% contrast
  condVar <- sigma - (cSig %*% t(cSig)) / sd^2
  sqrtCondVar <- sqrtMat(condVar)
  samp <- matrix(rnorm(nsamp * p), nrow = p, ncol = nsamp)
  samp <- sqrtCondVar %*% samp
  condMean <- mu + cSig %*% t(contSamp - m) / sd^2
  samp <- samp + condMean
  return(t(samp))
}
