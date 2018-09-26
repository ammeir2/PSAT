# # Cluster arguments
# args <- commandArgs(TRUE)
# eval(parse(text=args[[1]]))
setting <- as.numeric(setting)

# Packages
library(PSAT)
library(mvtnorm)
library(expm)
library(magrittr)
library(MASS)

# Helper functions
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

runSim <- function(config, seed) {
  set.seed(seed)
  n <- config[["n"]]
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
  fdr <- matrix(nrow = reps, ncol = 5)
  power <- matrix(nrow = reps, ncol = 5)
  tries <- 0
  for(m in 1:reps) {
    # Generating data
    X <- rnorm(n * p) %>% matrix(nrow = n)
    X <- X %*% sqrtmat
    X <- scale(X)
    nonzero <- sample.int(p, sparsity)
    true <- rep(0, p)
    true[nonzero] <- rnorm(sparsity)
    if(snr == 0) {
      true[nonzero] <- 0
    }
    mu <- as.numeric(X %*% true)
    if(snr == 0) {
      ysig <- 1
    } else {
      ysig <- sqrt(var(mu)) / (snr * sqrt(400 / n))
    }
    waldPval <- 1
    while(waldPval > pvalThreshold) {
      tries <- tries + 1
      y <- rnorm(n, mu, ysig)
      lmfit <- lm(y ~ X - 1)
      empCov <- vcov(lmfit)
      empInvCov <- solve(empCov)
      naive <- coef(lmfit)
      wald <- as.numeric(t(naive) %*% empInvCov %*% naive)
      waldPval <- 1 - pchisq(wald, df = p)
      # cat(waldPval, " ")
    }
    # cat("\n")
    # cat("selection prob: ", m * 1/tries, "\n")
    threshold <- qchisq(1 - pvalThreshold, df = p)
    naiveSD <- sqrt(sum(lmfit$residuals^2) / (n - p))
    if(waldPval > pvalThreshold * 0.05^2) {
      nullSD <- sd(y)
    } else {
      nullSD <- naiveSD
    }
    nullCov <- empCov * nullSD^2 / naiveSD^2
    knownCov <- empCov * ysig^2 / naiveSD^2
    # knownCov <- covar * ysig^2

    # PSAT Analysis
    ctrl <- psatControl(nSamples = 5000, truncPmethod = "symmetric", quadraticSampler = "tmg")
    time <- system.time(psatfit <- mvnQuadratic(naive, nullCov, testMat = empInvCov, threshold = threshold,
                            verbose = FALSE, control = ctrl))[3]
    # cat("PSAT time: ", time, "\n")
    
    # known variance
    knownVar <- mvnQuadratic(naive, knownCov, testMat = empInvCov, threshold = threshold,
                             verbose = FALSE, control = ctrl)
    
    # naive variance 
    # naiveVar <- mvnQuadratic(naive, empCov, testMat = "wald", threshold = threshold,
    #                          verbose = FALSE, control = ctrl)
    
    # p-values
    hybrid <- getPval(psatfit, type = "hybrid")
    polyhedral <- getPval(psatfit, type = "polyhedral")
    naive <- 2 * pnorm(-abs(naive / sqrt(diag(empCov))))
    knownHybrid <- getPval(knownVar, type = "hybrid")
    knownPolyhedral <- getPval(knownVar, type = "polyhedral")
    # naiveHybrid <- getPval(naiveVar, type = "hybrid")
    # naivePolyhedral <- getPval(naiveVar, type = "polyhedral")
    pvals <- cbind(hybrid = hybrid, polyhedral = polyhedral, naive = naive, 
                   knownHybrid = knownHybrid, knownPolyhedral = knownPolyhedral,
                   # naiveHybrid = naiveHybrid, naivePolyhedral = naivePolyhedral,
                   true = true)
    fdr[m, ] <- apply(pvals[, 1:5], 2, compFDR, true)
    power[m, ] <- apply(pvals[, 1:5], 2, compPower, true)
    
    # Reporting intermediate
    # print(c(m, config))
    # print(colMeans(fdr[1:m, , drop = FALSE]))
    # print(colMeans(power[1:m, , drop = FALSE]))
    result[[1 + m]] <- pvals
  }
  
  result[[1]][["R_squared"]] <- var(mu) / (var(mu) + ysig^2)
  # print(result[[1]])
  return(result)
}

configurations <- expand.grid(n = c(200, 100, 50),
                              p = c(20), 
                              snr = c(0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2),
                              # snr = 0,
                              sparsity = c(3), 
                              rho = c(0, 0.5),
                              seed = 1,
                              pvalThreshold = 0.01,
                              reps = 1)
# configurations <- configurations[sample(nrow(configurations)), ]
# configurations <- configurations[1:2, ]

# configurations <- expand.grid(n = c(100, 200), 
#                               p = c(10), 
#                               snr = c(0.01, 0.02), 
#                               sparsity = c(1, 3), 
#                               rho = c(0, 0.5),
#                               seed = 1,
#                               pvalThreshold = 0.01,
#                               reps = 2)

library(foreach)
library(doParallel)
registerDoParallel(cores = 1)
foreach(i = 601:800) %dopar% {
  print(i)
  system.time(results <- apply(configurations, 1, runSim, seed = i))
  filename <- paste("simulations/results/unknownVar_FIXED_B_seed_", i, ".rds", sep = "")
  saveRDS(results, file = filename)
}
# filename <- paste("results/unknownVar_wKnownNaive_A_seed", setting, ".rds", sep = "")
# saveRDS(results, file = filename)
# saveRDS(results, file = "simulations/results/unknownVar_zeroSignal_A2.rds")
