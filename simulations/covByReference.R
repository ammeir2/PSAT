# Cluster arguments
# args <- commandArgs(TRUE)
# eval(parse(text=args[[1]]))
# setting <- as.numeric(setting)

# Packages
library(PSAT)
library(mvtnorm)
library(expm)
library(magrittr)

# Helper functions
checkSingular <- function(x) {
  inv <- NULL
  try(inv <- solve(x))
  if(is.null(inv)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

genMAF <- function(nVariants, gammaShape, gammaRate, minMAF, maxMAF) {
  mafs <- numeric(nVariants)
  for(i in 1:length(mafs)) {
    samp <- Inf
    while(samp < minMAF | samp > maxMAF) {
      samp <- rgamma(1, gammaShape, gammaRate)
    }
    mafs[i] <- samp
  }
  return(mafs)
}

genX <- function(sqrtMat, nSubjects,
                 minMAF, MAFthreshold, 
                 gammaShape = 1, gammaRate = 50, seed = NULL) {
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  nVariants <- nrow(sqrtMat)
  n <- nSubjects
  normGene <- rnorm(2 * n * nVariants) %>% matrix(ncol = nVariants)
  normGene <- normGene %*% sqrtMat
  mafs <- genMAF(nVariants, gammaShape = gammaShape, 
                 gammaRate = gammaRate, 
                 minMAF = 2 / nSubjects, maxMAF = MAFthreshold)
  X <- apply(normGene, 1, function(x) pnorm(x) < mafs) %>% t()
  X <- X[1:nSubjects, ] + X[(nSubjects + 1):(2 * nSubjects), ]
  someNotZero <- apply(X, 2, function(x) any(x != 0))
  X <- X[, someNotZero]
  nVariants <- ncol(X)
  sparseX <- Matrix(X, sparse = TRUE)
  X <- scale(X)
  return(X)
}

compFDR <- function(pvals, true, level = 0.05) {
  qvals <- p.adjust(pvals, method = "BH")
  if(any(qvals < level)) {
    return(sum(true == 0 & qvals < level) / sum(qvals < level))
  } else {
    return(0)
  }
}

compPower <- function(pvals, true, level = 0.05) {
  qvals <- p.adjust(pvals, method = "BH")
  return(sum(true != 0 & qvals < level) / sum(true != 0))
}

runSim <- function(config, seed, verbose = TRUE) {
  set.seed(seed)
  n <- config[["n"]]
  refSize <- config[["refSize"]]
  p <- config[["p"]]
  snr <- config[["snr"]]
  sparsity <- config[["sparsity"]]
  rho <- config[["rho"]]
  pvalThreshold <- config[["pvalThreshold"]]
  reps <- config[["reps"]]
  
  covar <- rho^as.matrix(dist(1:p))
  invcov <- solve(covar)
  sqrtmat <- expm(covar)
  
  result <- vector(reps + 1, mode = "list")
  result[[1]] <- config
  fdr <- matrix(nrow = reps, ncol = 3)
  power <- matrix(nrow = reps, ncol = 3)
  tries <- 0
  for(m in 1:reps) {
    # Generating data
    refSingular <- TRUE
    while(refSingular) {
      X <- genX(sqrtmat, n + refSize, 
                minMAF = 2 / n, MAFthreshold = 0.1)
      refX <- X[1:refSize, ]
      refSingular <- checkSingular(var(refX))
    }
    X <- X[(refSize + 1):nrow(X), ]
    refX <- scale(refX)
    X <- scale(X)
    true <- rep(0, p)
    if(snr == 0) {
      ysig <- 1
      mu <- rep(0, n)
    } else {
      nonzero <- sample.int(p, sparsity)
      true[nonzero] <- rnorm(sparsity)
      mu <- as.numeric(X %*% true)
      ysig <- sqrt(var(mu)) / snr
    }
    waldPval <- 1
    refXtX <- t(refX) %*% refX
    scaledRefXtX <- refXtX * n / refSize
    refXtXinv <- solve(scaledRefXtX)
    while(waldPval > pvalThreshold) {
      tries <- tries + 1
      y <- rnorm(n, mu, ysig)
      y <- y - mean(y)
      suffStat <- t(X) %*% y
      varEst <- as.numeric(sum(y * y) - t(suffStat) %*% refXtXinv %*% suffStat) / n
      naive <- as.numeric(refXtXinv %*% t(X) %*% y)
      wald <- as.numeric(t(naive) %*% scaledRefXtX %*% naive) / varEst
      waldPval <- 1 - pchisq(wald, df = p)
    }
    if(verbose) cat("selection prob: ", m * 1/tries, "\n")
    threshold <- qchisq(1 - pvalThreshold, df = p)
    naiveSD <- sqrt(varEst)
    if(waldPval < pvalThreshold * 0.05^2) {
      nullSD <- naiveSD
    } else {
      nullSD <- sd(y)
    }
    
    refCov <- refXtXinv * nullSD^2
    refNaiveCov <- refXtXinv * varEst
    testMat <- scaledRefXtX / varEst
    
    # PSAT Analysis
    ctrl <- psatControl(nSamples = 10^4, quadraticSampler = "tmg", truncPmethod = c("symmetric"))
    threshold <- qchisq(1 - pvalThreshold, df = 50)
    time <- system.time(psatfit <- mvnQuadratic(naive, refCov, testMat = testMat,
                                                estimate = "naive",
                                                ci_type = c("naive"),
                                                threshold = threshold,
                                                verbose = FALSE, control = ctrl))[3]
    if(verbose) cat("PSAT time: ", time, "\n")
    
    # p-values
    hybrid <- getPval(psatfit, type = "hybrid")
    polyhedral <- getPval(psatfit, type = "polyhedral")
    naive <- 2 * pnorm(-abs(naive / sqrt(diag(refNaiveCov))))
    pvals <- cbind(hybrid = hybrid, polyhedral = polyhedral, naive = naive, true = true)
    fdr[m, ] <- apply(pvals[, 1:3], 2, compFDR, true)
    power[m, ] <- apply(pvals[, 1:3], 2, compPower, true)
    
    # Reporting intermediate
    if(verbose) {
      print(c(m, config))
      print(colMeans(fdr[1:m, , drop = FALSE]))
      print(colMeans(power[1:m, , drop = FALSE]))
    }
    result[[1 + m]] <- pvals
  }
  
  result[[1]][["R_squared"]] <- var(mu) / (var(mu) + ysig^2)
  if(verbose) print(result[[1]])
  return(result)
}

configA <- expand.grid(n = c(10^4),
                       refSize = c(625, 1250, 2500, 5000, 10^4),
                       p = c(50),
                       snr = c(0.015, 0.0315, 0.0625, 0.09, 0.125),
                       sparsity = c(1, 3, 6),
                       rho = c(0),
                       seed = 1,
                       pvalThreshold = 0.001,
                       reps = 1)
configB <- expand.grid(n = c(10^4),
                       refSize = c(1250, 2500, 5000, 10^4, 2 * 10^4),
                       p = c(50),
                       snr = c(0.125, 0.2, 0.25, 0.35, 0.5, 0.7),
                       sparsity = c(1, 3, 6),
                       rho = c(0.8),
                       seed = 1,
                       pvalThreshold = 0.001,
                       reps = 1)
configC <- expand.grid(n = c(10^4),
                       refSize = c(625, 1250, 2500, 5000, 10^4, 2 * 10^4),
                       p = c(50), 
                       snr = c(0.015, 0.0315, 0.0625, 0.09, 0.125, 0.2, 0.25, 0.35, 0.5, 0.7), 
                       sparsity = c(1, 3, 6), 
                       rho = c(0.4),
                       seed = 1,
                       pvalThreshold = 0.001,
                       reps = 1)
configurations <- rbind(configC, configA, configB)
setting <- setting

library(foreach)
library(doParallel)
registerDoParallel(cores = 1)
foreach(i = 201:300) %dopar% {
  print(i)
  system.time(results <- apply(configurations, 1, runSim, seed = i, verbose = FALSE))
  filename <- paste("simulations/results/covByRef_FIXED_A_seed_", i, ".rds", sep = "")
  saveRDS(results, file = filename)
}

# system.time(results <- apply(configurations, 1, runSim, seed = setting))
# filename <- paste("results/covByRef_H_seed", setting, ".rds", sep = "")
# saveRDS(results, file = filename)
