quadraticRB <- function(y, sigma, testmat, threshold,
                        alpha, rbIters = NULL,
                        variables = NULL,
                        computeFull = FALSE) {
  p <- length(y)
  if(is.null(variables)) {
    variables <- diag(length(y))
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

  ci <- matrix(nrow = nrow(variables), ncol = 2)
  for(j in 1:nrow(variables)) {
    contrast <- variables[j, ]
    contsd <- sqrt(as.numeric(t(contrast) %*% sigma %*% contrast))
    observed <- sum(contrast * y)
    
    for(upper in c(TRUE, FALSE)) {
      if(!computeFull) {
        if(upper & observed > 0) {
          next
        } else if(!upper & observed < 0) {
          next
        }
      }

      if(upper) {
        init <- observed + qnorm(1 - alpha / 4) * contsd
        a <- alpha / 2
      } else {
        init <- observed - qnorm(1 - alpha / 4) * contsd
        a <- alpha / 2
      }

      adjVecs <- computeAdjVecs(y, contrast, sigma)
      u <- quadraticRobinsMonroe(eta = contrast, 
                                 addVec = adjVecs$addVec, multVec = adjVecs$multVec,
                                 upper = upper, alpha = a,
                                 observed = observed, initU = init,
                                 threshold = threshold,
                                 sqrtTestMat = sqrtTestMat,
                                 invSqrtTestMat = invSqrtTestMat,
                                 sampPrecision = sampP,
                                 stepSize = 17 * sqrt(diag(sigma)[j]),
                                 burnin = 1000, mhiters = 5, rbIters = rbIters)
      ci[j, upper + 1] <- u
    }
  }

  ci[is.na(ci[, 1]), 1] <- -Inf
  ci[is.na(ci[, 2]), 2] <- Inf

  return(ci)
}

#' @import magrittr
linearRB <- function(y, sigma, testVec, threshold,
                     alpha, rbIters = NULL,
                     variables = NULL,
                     computeFull = FALSE) {
  p <- length(y)
  if(is.null(variables)) {
    variables <- diag(p)
  }

  # Computing matrices ---------------
  if(is.null(rbIters)) {
    rbIters <- 500  * (1 - alpha) / alpha
  }

  testsd <- as.numeric(sqrt(t(testVec) %*% sigma %*% testVec))
  ci <- matrix(nrow = nrow(variables), ncol = 2)
  for(j in 1:nrow(variables)) {
    contrast <- variables[j, ]
    observed <- sum(contrast * y)
    contsd <- as.numeric(sqrt(t(contrast) %*% sigma %*% contrast))
    
    for(upper in c(TRUE, FALSE)) {
      if(!computeFull) {
        if(upper & observed > 0) {
          next
        } else if(!upper & observed < 0) {
          next
        }
      }

      if(upper) {
        init <- observed + qnorm(1 - alpha / 4) * contsd
        a <- alpha / 2
      } else {
        init <- observed - qnorm(1 - alpha / 4) * contsd
        a <- alpha / 2
      }

      condSD <- contsd^2 - (t(contrast) %*% sigma %*% testVec)^2 / testsd^2 
      condSD <- as.numeric(sqrt(condSD))
      regConst <- as.numeric(t(testVec) %*% sigma %*% contrast) / testsd^2
      adjVecs <- computeAdjVecs(y, contrast, sigma) %>% sapply(function(x) sum(testVec * x))
      u <- linearRobinsMonroe(upper = upper, alpha = alpha, 
                              observed = observed, 
                              addVec = adjVecs[1], multVec = adjVecs[2],
                              initU = init,
                              threshold = threshold,
                              contsd = testsd, condSD = condSD, 
                              regConst = regConst,
                              stepSize = 17 * contsd,
                              rbIters = rbIters)
      ci[j, upper + 1] <- u
    }
  }

  ci[is.na(ci[, 1]), 1] <- -Inf
  ci[is.na(ci[, 2]), 2] <- Inf

  return(ci)
}


computeAdjVecs <- function(y, eta, sigma) {
  p <- length(y)
  proj <- diag(p) - eta %*% solve(t(eta) %*% eta) %*% t(eta)
  basis <- matrix(rnorm((p - 1) * p), nrow = p - 1)
  basis <- basis %*% proj
  basis <- rbind(eta, basis)
  
  regMat <- sigma %*% t(basis) %*% solve(basis %*% sigma %*% t(basis))
  addVec <- y - regMat %*% basis %*% y
  multVec <- regMat[, 1]
  return(list(addVec = addVec, multVec = multVec))
}

