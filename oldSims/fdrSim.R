args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)

library(ggplot2)
library(dplyr)
library(PSAT)
library(MASS)

run.simulation <- function(config) {
  gen_genotype = function(Maf, N, sqrtSig){
    J <- length(Maf)
    randA <- matrix(rnorm(N * length(Maf)), nrow = N) %*% sqrtSig
    randB <- matrix(rnorm(N * length(Maf)), nrow = N) %*% sqrtSig
    dat <- randA
    for(i in 1:J) {
      dat[, i] <- (randA[, i] <= quantile(randA[, i], Maf[i])) + (randB[, i] <= quantile(randB, Maf[i]))
    }

    return(dat)
  }

  n <- config[["n"]]
  p <- config[["p"]]
  signal <- config[["signal"]]
  sparsity <- config[["sparsity"]]
  reps <- config[["reps"]]
  betaType <- config[["betaType"]]
  noiseType <- config[["noiseType"]]
  pthreshold <- config[["pthreshold"]]
  MAFthreshold <- config[["MAFthreshold"]]
  rho <- config[["rho"]]
  ysig <- 1

  sigma <- rho^as.matrix(dist(1:p, method = "euclidean"))
  sqrtSig <- expm::sqrtm(sigma)

  attempt <- 0
  maxAttempts <- 10^3 * reps
  simResults <- list()
  allNaiveNull <- c()
  allMleNull <- c()
  allScoreNull <- c()
  allPolyNull <- c()
  for(m in 1:reps) {
    iterResults <- list()
    print(round(c(rep = m, n = n, p = p, signal = signal, sparsity = sparsity, pselect = (m - 1) / attempt), 3))
    print(c(threshold = pthreshold))
    hat <- NA
    while(is.na(hat[1])) {
      AAA <- rgamma(300, 1, 300)
      MAF <- (AAA[AAA>0.0002 & AAA < MAFthreshold])[1:p]
      X <- gen_genotype(MAF, n, sqrtSig)
      X <- apply(X, 2, function(x) x - mean(x))
      XtX <- t(X) %*% X
      if(betaType == 0) {
        trueBeta <- c(rep(signal, sparsity), rep(0, p - sparsity)) / sparsity
      } else if(betaType == 1) {
        trueBeta <- rexp(sparsity, 1)
        trueBeta <- trueBeta / sum(trueBeta) * signal
        trueBeta <- trueBeta * (1 - 2*rbinom(sparsity, 1, 0.5))
        trueBeta <- c(trueBeta, rep(0, p - sparsity))
        trueBeta <- trueBeta[order(runif(p))]
      }
      trueMu <- X %*% trueBeta
      trueBeta <- trueBeta / sd(trueMu) * signal
      trueBeta[is.nan(trueBeta)] <- 0
      trueMu <- X %*% trueBeta
      try(hat <- solve(t(X) %*% X) %*% t(X))
    }

    wPval <- 1
    threshold <- pthreshold
    chisqThreshold <- qchisq(threshold, df = ncol(X), lower.tail = FALSE)
    iterAttempt <- 0
    while(wPval > threshold) {
      iterAttempt <- iterAttempt + 1
      attempt <- attempt + 1
      if(noiseType == 0) {
        y <- rnorm(n, trueMu, sd = ysig)
      } else if(noiseType == 1) {
        df <- 7
        noise <- rt(n, 7) * sqrt((df - 2) / df)
        y <- trueMu + noise
      } else if(noiseType == 2) {
        noise <- rexp(n, sqrt(2)) * (1 - 2 * rbinom(n, 1, 0.5))
        y <- trueMu + noise
      } else if(noiseType == 3) {
        c <- sqrt(12)/2
        noise <- runif(n, min = -c, max = c)
        y <- trueMu + noise
      }
      y <- y - mean(y)
      mleBeta <- as.vector(hat %*% y)
      wald <- t(mleBeta) %*% XtX %*% mleBeta / ysig^2
      wPval <- pchisq(wald, ncol(X), lower.tail = FALSE)
      if(attempt > maxAttempts) break
      if(attempt %% 1000 == 0) print(paste("ATTEMPT = ", attempt))
      if(iterAttempt > 10^5) {
        iterAttempt <- 0
        while(is.na(hat[1])) {
          AAA = rgamma(300, 1, 300)
          MAF = (AAA[AAA>0.0002 & AAA<0.01])[1:p]
          X <- gen_genotype(MAF,n)
          X <- apply(X, 2, function(x) x - mean(x))
          XtX <- t(X) %*% X
          trueBeta <- c(rep(signal, sparsity), rep(0, p - sparsity)) / sparsity
          trueMu <- X %*% trueBeta
          try(hat <- solve(XtX) %*% t(X))
        }
      }
    }
    if(attempt > maxAttempts) break

    # Computing fit ------------------------
    ctrl <- psatControl(truncPmethod = "symmetric")
    allfit <- mvnQuadratic(y = mleBeta, sigma = solve(XtX), threshold = chisqThreshold,
                           pvalue_type = c("hybrid", "polyhedral", "naive", "global_null"),
                           ci_type = c("naive"),
                           verbose = FALSE, control = ctrl)

    zero <- trueBeta == 0
    nonzero <- trueBeta != 0

    # Saving results -------------------------
    iterResults$params <- c(m = m, n = n, p = p, sparsity = sparsity,
                            signal = signal,
                            selectProb = 1 / iterAttempt, chisqPval = wPval,
                            coefType = betaType, noiseType = noiseType,
                            MAFthreshold = MAFthreshold, pthreshold = pthreshold,
                            R_squared = var(trueMu) / (1 + var(trueMu)))
    iterResults$naiveNullPvals <- getPval(allfit, type = "naive")[zero]
    iterResults$naiveAltPvals <- getPval(allfit, type = "naive")[nonzero]
    iterResults$conditionalNullPvals <- getPval(allfit, type = "global-null")[zero]
    iterResults$conditionalAltPvals <- getPval(allfit, type = "global-null")[zero]
    iterResults$nullNotScore <- getPval(allfit, type = "hybrid")[zero]
    iterResults$altNotScore <- getPval(allfit, type = "hybrid")[nonzero]
    iterResults$polyNull <- getPval(allfit, type = "polyhedral")[zero]
    iterResults$polyAlt <- getPval(allfit, type = "polyhedral")[nonzero]

    allNaiveNull <- c(allNaiveNull, iterResults$naiveNullPvals)
    allMleNull <- c(allMleNull, iterResults$conditionalNullPvals)
    allScoreNull <- c(allScoreNull, iterResults$nullNotScore)
    allPolyNull <- c(allPolyNull, iterResults$polyNull)

    simResults[[m]] <- iterResults
  }

  return(simResults)
}

#FOR FDR ----------------------------
configurations <- expand.grid(n = c(10^4),
                              p = c(50),
                              noiseType = c(0, 2, 3),
                              signal = c(10^-3 * 2^(5:6), 0),
                              sparsity = c(3),
                              betaType = 1,
                              pthreshold = 10^-3,
                              MAFthreshold = 0.1,
                              reps = 2,
                              rho = 0.8)

set.seed(setting)
system.time(simResults <- apply(configurations, 1, run.simulation))
filename <- paste("results/fdr_sim_A_seed", setting, ".rds", sep = "")
saveRDS(simResults, file = filename)

