getQudraticLam <- function(testMat, sigma) {
  c <- chol(sigma)
  tempmat <- c %*% testMat %*% t(c)
  eig <- eigen(tempmat)
  vec <- eig$vectors
  P <- t(vec)
  lam <- eig$values
  return(lam)
}

getQuadraticThreshold <- function(quantile, lam) {
  lower <- 0
  upper <- 1
  prob <- CompQuadForm::liu(upper, lam)
  while(prob > quantile) {
    lower <- upper
    upper <- upper * 2
    prob <- CompQuadForm::liu(upper, lam)
  }

  q <- uniroot(function(x) CompQuadForm::liu(x, lam) - quantile,
               interval = c(lower, upper))$root
  return(q)
}

getTestThreshold <- function(alpha, sigma, testMat) {
  return(getQuadraticThreshold(alpha, getQudraticLam(testMat, sigma)))
}
