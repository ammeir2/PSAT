quadraticRB <- function(y, sigma, testmat, threshold,
                        alpha, rbIters = NULL,
                        variables = NULL,
                        computeFull = FALSE) {
  p <- length(y)
  if(is.null(variables)) {
    variables <- 1:p
  }

  # Computing matrices ---------------
  testEigen <- eigen(testmat)
  invTestMat <- testEigen$vectors %*% diag(1/(testEigen$values))
  invTestMat <- invTestMat %*% t(testEigen$vectors)
  sqrtTestMat <- testEigen$vectors %*% diag(sqrt(testEigen$values))
  sqrtTestMat <- sqrtTestMat %*% t(testEigen$vectors)
  invSqrtTestMat <- testEigen$vectors %*% diag(1 / sqrt(testEigen$values))
  invSqrtTestMat <- invSqrtTestMat %*% t(testEigen$vectors)
  sampSig <- t(sqrtTestMat) %*% sigma %*% sqrtTestMat
  sampP <- solve(sampSig)

  if(is.null(rbIters)) {
    rbIters <- (50 * (1 - alpha) / alpha)
  }

  ci <- matrix(nrow = p, ncol = 2)
  for(j in variables) {
    for(upper in c(TRUE, FALSE)) {
      if(!computeFull) {
        if(upper & y[j] > 0) {
          next
        } else if(!upper & y[j] < 0) {
          next
        }
      }

      if(upper) {
        init <- y[j] + qnorm(1 - alpha / 4) * sqrt(diag(sigma)[j])
        a <- alpha / 2
      } else {
        init <- y[j] - qnorm(1 - alpha / 4) * sqrt(diag(sigma)[j])
        a <- alpha / 2
      }

      u <- quadraticRobinsMonroe(var = j, upper = upper, alpha = a,
                                 observed = y[j], initU = init,
                                 threshold = threshold,
                                 sqrtTestMat = sqrtTestMat,
                                 invSqrtTestMat = invSqrtTestMat,
                                 sampPrecision = sampP,
                                 stepSize = 17 * sqrt(diag(sigma)[j]),
                                 burnin = 1000, mhiters = 5, rbIters = rbIters)
      ci[j, upper + 1] <- u
    }
    # print(ci[j, ])
  }

  ci[is.na(ci[, 1]), 1] <- -Inf
  ci[is.na(ci[, 2]), 2] <- Inf

  return(ci)
}


linearRB <- function(y, sigma, contrast, threshold,
                        alpha, rbIters = NULL,
                        variables = NULL,
                        computeFull = FALSE) {
  p <- length(y)
  if(is.null(variables)) {
    variables <- 1:p
  }

  # Computing matrices ---------------
  if(is.null(rbIters)) {
    rbIters <- 500  * (1 - alpha) / alpha
  }

  contsd <- as.numeric(sqrt(t(contrast) %*% sigma %*% contrast))
  ci <- matrix(nrow = p, ncol = 2)
  for(j in variables) {
    for(upper in c(TRUE, FALSE)) {
      if(!computeFull) {
        if(upper & y[j] > 0) {
          next
        } else if(!upper & y[j] < 0) {
          next
        }
      }

      if(upper) {
        init <- y[j] + qnorm(1 - alpha / 4) * sqrt(diag(sigma)[j])
        a <- alpha / 2
      } else {
        init <- y[j] - qnorm(1 - alpha / 4) * sqrt(diag(sigma)[j])
        a <- alpha / 2
      }

      condSD <- as.numeric(sqrt(sigma[j, j] - sum(sigma[j, ] * contrast)^2 / contsd^2))
      regConst <- sum(sigma[j, ] * contrast) / contsd^2
      u <- linearRobinsMonroe(upper, alpha, y[j], init,
                              contrast[j], threshold,
                              contsd, condSD, regConst,
                              stepSize = 17 * sqrt(diag(sigma)[j]),
                              rbIters = rbIters)
      ci[j, upper + 1] <- u
    }
  }

  ci[is.na(ci[, 1]), 1] <- -Inf
  ci[is.na(ci[, 2]), 2] <- Inf

  return(ci)
}


