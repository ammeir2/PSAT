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
  rho <- config[["rho"]]
  ysig <- 1

  sigma <- rho^as.matrix(dist(1:p, method = "euclidean"))
  sqrtSig <- expm::sqrtm(sigma)

  pselect <- 0
  simResults <- list()
  allNaiveNull <- c()
  allMleNull <- c()
  allScoreNull <- c()
  allPolyNull <- c()
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
      trueMu <- X %*% trueBeta
      trueBeta <- trueBeta / sd(trueMu) * signal
      trueBeta[is.nan(trueBeta)] <- 0
      trueMu <- X %*% trueBeta
      try(betaVar <- solve(XtX))
    }

    wPval <- 1
    threshold <- pthreshold
    chisqThreshold <- qchisq(threshold, df = ncol(X), lower.tail = FALSE)
    ncp <- t(trueBeta) %*% XtX %*% trueBeta / ysig^2
    pselect <- pchisq(chisqThreshold, df = ncol(X), ncp = ncp, lower.tail = FALSE)

    betaVar <- betaVar * ysig^2
    samplingBeta <- trueBeta
    naiveBeta <- sampleQuadraticConstraint(trueBeta, betaVar,
                                           chisqThreshold, XtX,
                                           sampSize = 1,
                                           burnin = 10^4)

    ctrl <- psatControl(truncPmethod = "symmetric")
    fit <- mvnQuadratic(naiveBeta, betaVar, threshold = chisqThreshold,
                        pvalue_type = c("naive", "polyhedral", "hybrid", "global-null"),
                        ci_type = c("naive"),
                        verbose = FALSE, control = ctrl)


    # Saving results -------------------------
    iterResults$params <- c(m = m, n = n, p = p, sparsity = sparsity,
                            signal = signal,
                            selectProb = pselect, chisqPval = wPval,
                            coefType = betaType, noiseType = noiseType,
                            MAFthreshold = MAFthreshold, pthreshold = pthreshold,
                            R_squared = var(trueMu) / (1 + var(trueMu)))
    zero <- trueBeta == 0
    nonzero <- trueBeta != 0
    iterResults$naiveNullPvals <- getPval(fit, type = "naive")[zero]
    iterResults$naiveAltPvals <- getPval(fit, type = "naive")[nonzero]
    iterResults$conditionalNullPvals <- getPval(fit, type = "global-null")[zero]
    iterResults$conditionalAltPvals <- getPval(fit, type = "global-null")[nonzero]
    iterResults$nullNotScore <- getPval(fit, type = "hybrid")[zero]
    iterResults$altNotScore <- getPval(fit, type = "hybrid")[nonzero]
    iterResults$polyNull <- getPval(fit, type = "polyhedral")[zero]
    iterResults$polyAlt <- getPval(fit, type = "polyhedral")[nonzero]

    allNaiveNull <- c(allNaiveNull, iterResults$naiveNullPvals)
    allMleNull <- c(allMleNull, iterResults$conditionalNullPvals)
    allScoreNull <- c(allScoreNull, iterResults$nullNotScore)
    allPolyNull <- c(allPolyNull, iterResults$polyNull)

    simResults[[m]] <- iterResults
  }

  return(simResults)
}

configurations <- expand.grid(n = c(10^4), p = 50,
                              signal = 2^(8:2) * 0.001,
                              sparsity = c(1, 2, 4, 8),
                              reps = 2,
                              betaType = 1,
                              noiseType = 0,
                              pthreshold = c(10^-3),
                              MAFthreshold = c(0.1),
                              rho = 0.8)

set.seed(setting)
system.time(simResults <- apply(configurations, 1, run.simulation))
filename <- paste("results/power_sim_A_seed", setting, ".rds", sep = "")
saveRDS(simResults, file = filename)



