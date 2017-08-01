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

#' Get the Threshold for a General Quadratic Test
#'
#' @description Computes the upper quantile of the distribution of a quadratic
#' form: \deqn{y' K y,} where \eqn{y ~ N(0, \Sigma)} and \eqn{K} is
#' positive semi-definite. This quantile is approximated using the
#' \code{\link[CompQuadForm]{liu}} method.
#'
#' @param alpha the upper quantile of the distribution/significance
#' level of the quadratic test.
#'
#' @param sigma the covariance of \eqn{y}.
#'
#' @param testMat the test matrix \eqn{K}.
#'
#' @return a numeric threshold for the quadratic test.
#'
#' @seealso \code{\link[CompQuadForm]{CompQuadForm}}, \code{\link[CompQuadForm]{liu}}
getTestThreshold <- function(alpha, sigma, testMat) {
  return(getQuadraticThreshold(alpha, getQudraticLam(testMat, sigma)))
}
