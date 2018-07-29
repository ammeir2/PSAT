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
    time <- system.time(condFit <- mvnQuadratic(naiveBeta, sigma, testMat = "wald",
                            pval_threshold = pthreshold,
                            pvalue_type = "naive",
                            ci_type = c("naive"),
                            confidence_level = 0.95, verbose = FALSE))
    print(time[3])


    # Coverage Rate ----------------------
    mle <- condFit$muhat
    naive <- condFit$y
    true <- trueBeta
    mle <- mean((mle - true)^2)
    naive <- mean((naive - true)^2)

    # Saving results -------------------------
    iterResults$params <- c(m = m, n = n, p = p, sparsity = sparsity,
                            signal = signal,
                            selectProb = pselect, chisqPval = 1,
                            coefType = betaType, noiseType = noiseType,
                            MAFthreshold = MAFthreshold, pthreshold = pthreshold,
                            t2type = t2type, mle = mle, naive = naive,
                            R_squared = var(trueMu) / (1 + var(trueMu)))

    simResults[[m]] <- iterResults
  }

  return(simResults)
}

configurations <- expand.grid(n = c(5000, 10^4, 1.5 * 10^4, 2 * 10^4),
                              p = c(5, 10, 20),
                              signal = c(0, 0.025, 0.05),
                              sparsity = c(2),
                              reps = 20,
                              betaType = 1,
                              noiseType = 0,
                              pthreshold = c(10^-3),
                              MAFthreshold = c(0.1),
                              t2type = c(2) ,
                              randRho = 0.8)

set.seed(setting)
system.time(simResults <- apply(configurations, 1, run.simulation))
filename <- paste("results/mle_sim_A_seed", setting, ".rds", sep = "")
saveRDS(simResults, file = filename)



