args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)

library(ggplot2)
library(dplyr)
library(PSAT)
library(MASS)

expit <- function(x) {
  1 / (1 + exp(-x))
}

logit <- function(x) {
  log(x / (1 - x))
}

run.simulation <- function(config) {
  gen_genotype = function(Maf, N, sqrtSig){
    J <- length(Maf)
    eta <- logit(Maf)
    rand <- matrix(rnorm(2 * N * length(Maf)), nrow = 2 * N) %*% sqrtSig
    dat <- rand
    for(i in 1:J) {
      dat[, i] <- rand[, i] <= quantile(rand[, i], Maf[i])
    }

    dat <- dat[1:N, ] + dat[(N+1):(2 * N), ]
    return(dat)
  }

  n <- config[[1]]
  p <- config[[2]]
  signal <- config[[3]]
  sparsity <- config[[4]]
  reps <- config[[5]]
  betaType <- config[[6]]
  noiseType <- config[[7]]
  pthreshold <- config[[8]]
  MAFthreshold <- config[[9]]
  t2type <- config[[10]]
  rho <- config[[11]]
  ysig <- 1

  sigma <- rho^as.matrix(dist(1:p, method = "euclidean"))
  sqrtSig <- expm::sqrtm(sigma)

  pselect <- 0
  simResults <- list()
  allNaiveNull <- c()
  allMleNull <- c()
  allScoreNull <- c()
  allPolyNull <- c()
  meancover <- 0
  meandet <- 0
  for(m in 1:reps) {
    iterResults <- list()
    print(round(c(rep = m, n = n, p = p, signal = signal, sparsity = sparsity, pselect = pselect), 3))
    print(c(threshold = pthreshold))
    betaVar <- NA
    while(is.na(betaVar[1])) {
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
      try(betaVar <- solve(XtX))
    }

    wPval <- 1
    threshold <- pthreshold
    chisqThreshold <- qchisq(threshold, df = ncol(X), lower.tail = FALSE)
    ncp <- t(trueBeta) %*% XtX %*% trueBeta / ysig^2
    pselect <- pchisq(chisqThreshold, df = ncol(X), ncp = ncp, lower.tail = FALSE)

    betaVar <- betaVar * ysig^2
    samplingBeta <- trueBeta
    naiveBeta <- sampleQuadraticConstraint(trueBeta, solve(XtX),
                                           chisqThreshold, XtX,
                                           sampSize = 1,
                                           burnin = 10^4)
    Xy <- as.numeric(XtX %*% t(naiveBeta))

    # Inference -------------------
    sigma <- solve(XtX)
    ctrl <- psatControl(truncPmethod = "symmetric")
    time <- system.time(condFit <- mvnQuadratic(naiveBeta, sigma, testMat = "wald",
                            pval_threshold = pthreshold,
                            pvalue_type = "polyhedral",
                            ci_type = c("polyhedral", "naive", "switch"),
                            confidence_level = 0.95, verbose = FALSE, control = ctrl))
    print(time[3])


    # Coverage Rate ----------------------
    tnCI <- getCI(condFit, type = "polyhedral")
    naiveCI <- getCI(condFit, type = "naive")
    condCI <- condFit$switchCI
    naiveCover <- mean((naiveCI[, 1] < trueBeta & naiveCI[, 2] > trueBeta))
    condCover <- mean((condCI[, 1] < trueBeta & condCI[, 2] > trueBeta))
    TNcover <- mean((tnCI[, 1] < trueBeta & tnCI[, 2] > trueBeta))

    # Power to determine sign ---------------
    zero <- trueBeta == 0
    nonzero <- trueBeta != 0
    naiveDet <- mean((naiveCI[nonzero, 1] > 0 | naiveCI[nonzero, 2] < 0))
    condDet <- mean((condCI[nonzero, 1] > 0 | condCI[nonzero, 2] < 0))
    TNdet <- mean((tnCI[nonzero, 1] > 0 | tnCI[nonzero, 2] < 0))

    # CI width -------------
    naiveSize <- mean(naiveCI[, 2] - naiveCI[, 1])
    condSize <- mean(condCI[, 2] - condCI[, 1])
    polySize <- mean(tnCI[, 2] - tnCI[, 1])

    # Saving results -------------------------
    iterResults$params <- c(m = m, n = n, p = p, sparsity = sparsity,
                            signal = signal,
                            selectProb = pselect, chisqPval = 1,
                            coefType = betaType, noiseType = noiseType,
                            MAFthreshold = MAFthreshold, pthreshold = pthreshold,
                            t2type = t2type, 
                            R_squared = var(trueMu) / (ysig^2 + var(trueMu)))
    iterResults$cover <- c(naive = naiveCover, cond = condCover, TN = TNcover)
    iterResults$det <- c(naive = naiveDet, cond = condDet, TN = TNdet)
    iterResults$naiveCI <- naiveCI
    iterResults$condCI <- condCI
    iterResults$tnCI <- tnCI
    iterResults$true <- trueBeta
    iterResults$realNaive <- condFit$naiveCI

    meancover <- meancover * (m - 1) / m + 1 / m * c(naiveCover, condCover, TNcover)
    meandet <- meandet * (m - 1) / m + 1 / m * c(naiveDet, condDet, TNdet)
    cat("cover ", meancover, "\n")
    cat("det ", meandet, "\n")

    simResults[[m]] <- iterResults
  }

  return(simResults)
}

configurations <- expand.grid(n = c(10^4), p = c(20),
                              signal = 2^(8:0) * 0.001,
                              sparsity = c(1, 4, 8, 2),
                              reps = 10,
                              betaType = 1,
                              noiseType = 0,
                              pthreshold = c(10^-3),
                              MAFthreshold = c(0.1),
                              t2type = c(2) ,
                              randRho = 0.8)
set.seed(setting)
system.time(simResults <- apply(configurations, 1, run.simulation))
filename <- paste("results/ci_sim_A_seed", setting, ".rds", sep = "")
saveRDS(simResults, file = filename)
