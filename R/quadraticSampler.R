sampleQuadraticConstraint <- function(mu, sigma, threshold, testMat,
                                      init = NULL,
                                      sampSize = 10^3,
                                      burnin = 1000,
                                      trim = 50) {
  if(threshold < 0) stop("Threshold must be larger than zero!")

  # Pre-processing
  testEigen <- eigen(testMat)
  invTestMat <- testEigen$vectors %*% diag(1/(testEigen$values))
  invTestMat <- invTestMat %*% t(testEigen$vectors)
  sqrtTestMat <- testEigen$vectors %*% diag(sqrt(testEigen$values))
  sqrtTestMat <- sqrtTestMat %*% t(testEigen$vectors)
  invSqrtTestMat <- testEigen$vectors %*% diag(1 / sqrt(testEigen$values))
  invSqrtTestMat <- invSqrtTestMat %*% t(testEigen$vectors)
  sigma <- t(sqrtTestMat) %*% sigma %*% sqrtTestMat
  precision <- solve(sigma)
  mu <- sqrtTestMat %*% mu
  if(!is.null(init)) {
    init <- sqrtTestMat %*% init
  } else {
    init <- mu
  }

  # Running C code
  sample <- sampleQuadraticMVTcpp(init, mu, precision, threshold,
                                  sampSize, burnin, trim)
  sample <- sample %*% invSqrtTestMat
  return(sample)
}
