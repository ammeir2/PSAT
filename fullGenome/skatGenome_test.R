# Information from cluster ---
# args <- commandArgs(TRUE)
# eval(parse(text=args[[1]]))
setting <- as.numeric(setting)

# Helper functions ----
library(magrittr)
library(PSAT)
library(MASS)
library(Matrix)
library(expm)

checkError <- function(pvals, nulls, method = "BH", level = 0.05) {
  qvals <- p.adjust(pvals, method = method)
  if(any(qvals < level)) {
    if(method == "BH") {
      error <- sum(nulls & qvals < level) / sum(qvals < level)
    } else if(method == "bonferroni") {
      error <- any(nulls & qvals < level)
    }
  }
  return(error)
}

checkInclusion <- function(ci, true) {
  inside <- 0
  for(i in 1:nrow(ci)) {
    if(ci[i, 1] < true[i] & ci[i, 2] > true[i]) {
      inside <- inside + 1
    }
  }
  
  return(inside / nrow(ci))
}

checkPower <- function(pvals, nulls, method = "BH", level = 0.05) {
  if(all(nulls)) {
    return(0)
  }
  
  qvals <- p.adjust(pvals, method = method)
  power <- sum(qvals < level & !nulls) / sum(!nulls)
  return(power)
}

# Data generating functions ------
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

genGene <- function(sqrtMats, nSubjects,
                    minMAF, MAFthreshold, notNull, 
                    sparsity,
                    snrLevels, snrProbs,
                    gammaShape = 1, gammaRate = 50, 
                    geneInd, seed = NULL) {
  if(!is.null(geneInd)) {
    set.seed(geneInd)
  }
  
  sqrtMat <- sqrtMats[[sample.int(length(sqrtMats), 1)]]
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
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  nonNull <- notNull[geneInd]
  trueCoef <- rep(0, nVariants)
  if(nonNull) {
    nonzero <- max(floor(nVariants * sparsity), 1)
    snr <- sample(snrLevels, 1, FALSE, snrProbs)
    probs <- 1 / mafs
    probs <- probs / sum(probs)
    whichNonZero <- sample.int(nVariants, nonzero, prob = probs)
    trueCoef[sample.int(nVariants, nonzero)] <- rnorm(nonzero)
    trueCoef <- trueCoef / sum(abs(trueCoef))
    mu <- as.numeric(X %*% trueCoef)
    mu <- mu / sd(mu) * snr
  } else {
    mu <- rep(0, nSubjects)
    nonzero <- 0
  }
  
  return(list(X = X, sparseX = sparseX, mu = mu, coef = trueCoef, nonzero = nonzero))
}

# Simulation function --------
runSim <- function(config, verbose = TRUE) {
  attach(config, warn.conflicts = FALSE)
  sqrtMats <- lapply(unlist(geneSizes), function(size) {
    covmat <- rho^as.matrix(dist(1:size))
    sqrtmat <- sqrtm(covmat)
    return(sqrtmat)
  })
  
  # Which genes are non Null?? ---
  set.seed(1)
  nNotNull <- ceiling(genomSize * (1 - pNullGenes))
  notNull <- c(rep(TRUE, nNotNull), rep(FALSE, genomSize - nNotNull)) %>% sample()
  
  # Seeds for genes -------
  set.seed(seed)
  seeds <- sample.int(10^6, genomSize)
  
  # SNR levels -----
  snrLevels <- c(1) 
  if(signal == "weak") {
    snrProbs <- 1
  } else if(signal == "strong") {
    snrProbs <- 1:5 / 15
  }
  
  # computing expected value ------
  yMean <- rep(0, nSubjects)
  print("Generating response")
  sparseMats <- list(genomSize, mode = "list")
  coefList <- list(genomSize, mode = "list")
  if(verbose) pb <- txtProgressBar(min = 0, max = genomSize, style = 3)
  totVariants <- 0  
  for(g in 1:genomSize) {
    # Generating data
    dat <- genGene(sqrtMats, nSubjects,
                   minMAF, MAFthreshold, notNull, 
                   sparsity,
                   snrLevels, snrProbs,
                   gammaShape = 1, gammaRate = 50,
                   geneInd = g, seed = seeds[g])  
    sparseMats[[g]] <- dat$sparseX
    coefList[[g]] <- dat$coef
    yMean <- yMean + dat$mu
    totVariants <- totVariants + ncol(dat$X)
    if(verbose) setTxtProgressBar(pb, g)
  }
  rm(dat)
  if(verbose) close(pb)
  yMean <- yMean / sqrt(var(yMean) * (1 - R) / R)
  y <- rnorm(nSubjects, mean = yMean, sd = 1)
  yMean <- yMean - mean(y)
  y <- y - mean(y)
  yvar <- var(y)
  config[["R_squared"]] <- var(yMean) / (var(yMean) + 1)
  print(var(yMean) / (var(yMean) + 1))
  
  # Aggregate Testing loop
  print("Aggregate Testing")
  slot <- 1
  set.seed(seed)
  selectedGenes <- vector(genomSize, mode = "list")
  if(verbose) pb <- txtProgressBar(min = 0, max = genomSize, style = 3)
  univPvals <- numeric(totVariants)
  trueCoefs <- numeric(totVariants)
  univInd <- 1
  for(g in 1:genomSize) {
    X <- as.matrix(sparseMats[[g]])
    trueCoef <- coefList[[g]]
    sparseMats[[g]] <- NA
    coefList[[g]] <- NA
    lmPvals <- summary(lm(y ~ X))$coefficients[-1, 4]
    for(i in 1:ncol(X)) {
      univPvals[univInd] <- lmPvals[i]
      trueCoefs[univInd] <- trueCoef[i]
      univInd <- univInd + 1
    }
    
    # Aggregate testing -------
    XtX <- t(X) %*% X
    XtXinv <- ginv(XtX)
    suffStat <- t(X) %*% y
    naiveCoef <- XtXinv %*% suffStat
    coefCov <- yvar * XtXinv
    waldMat <- XtX / as.numeric(var(y - X %*% naiveCoef))
    waldStat <- as.numeric(t(naiveCoef) %*% waldMat %*% naiveCoef)
    waldPval <- pchisq(waldStat, df = ncol(X), lower.tail = FALSE)
    critVal <- qchisq(pvalThreshold, df = ncol(X), lower.tail = FALSE)
    naiveSD <- sqrt(diag(coefCov)) * sd(y - as.numeric(X %*% naiveCoef)) / sd(y)
    trueProj <- XtXinv %*% t(X) %*% yMean
    # print(c(g, ncol(X), waldStat, critVal, sum(abs(trueCoef)), dat$nonzero))
    if(waldStat > critVal) {
      selectedGenes[[slot]] <- list(obs = naiveCoef, 
                                    cov = coefCov, 
                                    naiveSD = naiveSD,
                                    true = trueCoef,
                                    trueProj = trueProj,
                                    waldMat = waldMat,
                                    threshold = critVal,
                                    waldPval = waldPval,
                                    test = "wald")
      slot <- slot + 1
    }
    
    # SKAT Testing ----------
    props <- apply(X, 2, function(x) mean(x > 0))
    skatWeights <- diag(dbeta(props, 1, 25)^2)
    skatMat <- XtX %*% skatWeights %*% XtX
    skatMat <- cov2cor(skatMat)
    quadlam <- PSAT:::getQudraticLam(skatMat, coefCov)
    skatThreshold <- PSAT:::getQuadraticThreshold(pvalThreshold, quadlam)
    skatStat <- as.numeric(t(naiveCoef) %*% skatMat %*% naiveCoef)
    skatPval <- CompQuadForm::liu(skatStat, quadlam)
    print(c(g, waldPval, skatPval))
    if(skatPval < pvalThreshold) {
      selectedGenes[[slot]] <- list(obs = naiveCoef, 
                                    cov = coefCov, 
                                    naiveSD = naiveSD,
                                    true = trueCoef,
                                    trueProj = trueProj,
                                    waldMat = skatMat,
                                    threshold = skatThreshold,
                                    waldPval = skatPval,
                                    test = "skat")
      slot <- slot + 1
    }
    
    if(verbose) setTxtProgressBar(pb, g)
  }
  
  
  rm(X)
  if(verbose) close(pb)
  
  # selection with BH
  selectedGenes <- selectedGenes[1:(slot - 1)]
  isWald <- sapply(selectedGenes, function(x) x$test == "wald")
  isSkat <- sapply(selectedGenes, function(x) x$test == "skat")
  waldGenes <- list()
  skatGenes <- list()
  if(sum(isWald) > 0) {
    aggPvals <- sapply(selectedGenes[isWald], function(x) x$waldPval)
    aggQvals <- p.adjust(aggPvals, method = "BH", n = genomSize)
    nWaldSelected <- sum(aggQvals < pvalThreshold)
    waldGenes <- selectedGenes[isWald][aggQvals < pvalThreshold]
  }
  
  if(sum(isSkat) > 0) {
    aggPvals <- sapply(selectedGenes[isSkat], function(x) x$waldPval)
    aggQvals <- p.adjust(aggPvals, method = "BH", n = genomSize)
    nSkatSelected <- sum(aggQvals < pvalThreshold)
    skatGenes <- selectedGenes[isSkat][aggQvals < pvalThreshold]
  }
  bhFracWald <- length(waldGenes) / genomSize
  bhFracSkat <- length(skatGenes) / genomSize
  waldPvalThreshold <- bhFracWald * pvalThreshold
  skatPvalThreshold <- bhFracSkat * pvalThreshold
  selectedGenes <- c(waldGenes, skatGenes)

  # Conducting Inference
  inferenceResults <- vector(length(selectedGenes) + 1, mode = "list")
  if(length(selectedGenes) > 0) {
    totVariants <- sapply(selectedGenes, function(x) ncol(x$cov)) %>% sum()
  } else {
    totVariants <- 0
  }
  if(length(selectedGenes) > 0) {
    print("Conducting Inference")
    if(verbose) pb <- txtProgressBar(min = 0, max = length(selectedGenes), style = 3)
    for(g in 1:length(selectedGenes)) {
      gene <- selectedGenes[[g]]
      naive <- gene$obs %>% as.numeric()
      true <- gene$true %>% as.numeric()
      trueProj <- gene$trueProj
      sds <- sqrt(diag(gene$cov))
      naiveSD <- gene$naiveSD
      truePvals <- 2 * pnorm(-abs(true /sds))
      empNonzero <- sum(truePvals < 0.05)
      test <- gene[["test"]]
      
      if(test == "wald") {
        pthreshold <- waldPvalThreshold
      } else {
        pthreshold <- skatPvalThreshold
      }
      
      psatFit <- NULL
      control <- psatControl(nSamples = 10^4, quadraticSampler = "tmg")
      try(psatFit <- mvnQuadratic(naive, gene$cov, testMat = gene$waldMat, 
                                  pval_threshold = pthreshold, 
                                  contrasts = diag(length(naive)),
                                  estimate_type = "naive",
                                  pvalue_type = c("hybrid", "polyhedral"), 
                                  ci_type = c("switch", "polyhedral"),
                                  verbose = FALSE,
                                  control = control))
      if(is.null(psatFit)) {
        control <- psatControl(nSamples = 10^4, quadraticSampler = "tmg", truncPmethod = "symmetric")
        psatFit <- mvnQuadratic(naive, gene$cov, testMat = gene$waldMat, 
                                pval_threshold = pthreshold, 
                                contrasts = diag(length(naive)),
                                estimate_type = "naive",
                                pvalue_type = c("hybrid", "polyhedral"), 
                                ci_type = c("switch", "polyhedral"),
                                verbose = FALSE,
                                control = control)
      }
      
      # Cover Rate
      switchCover <- checkInclusion(psatFit$switchCI, trueProj)
      polyCover <- checkInclusion(psatFit$polyCI, trueProj)
      bbQuantile <- pthreshold / 2
      bbCritVal <- qnorm(1 - bbQuantile)
      bbCI <- naive + cbind(-naiveSD * bbCritVal, naiveSD * bbCritVal)
      bbCover <- checkInclusion(bbCI, trueProj)
      
      # FDR
      fdrFrame <- data.frame(hybrid = getPval(psatFit, type = "hybrid"), 
                             poly = getPval(psatFit, type = "polyhedral"), 
                             naive = getPval(psatFit, type = "naive"), 
                             bbFraction = pthreshold * pvalThreshold, 
                             true = true)
      # Avg Power 
      nulls <- true == 0
      bbLevel <- length(selectedGenes) / genomSize * 0.05
      hybridPower <- getPval(psatFit, type = "hybrid") %>% 
        checkPower(nulls, method = "BH", level = 0.05)
      polyPower <- getPval(psatFit, type = "polyhedral") %>% 
        checkPower(nulls, method = "BH", level = 0.05)
      bbPower <- getPval(psatFit, type = "naive") %>% 
        checkPower(nulls, method = "BH", level = bbLevel)
      
      # Reporting
      geneResult <- list()
      config[["test"]] <- test
      geneResult$config <- config
      geneResult$cover <- c(switch = switchCover, poly = polyCover, bb = bbCover)
      geneResult$power <- c(hybrid = hybridPower, poly = polyPower, bbAvgPower = bbPower)
      geneResult$FDR <- fdrFrame
      inferenceResults[[g]] <- geneResult
      
      if(verbose) setTxtProgressBar(pb, g)
    }
    if(verbose) close(pb)
  }
  
  if(length(inferenceResults) == 1) {
    inferenceResults[[1]] <- list()
    inferenceResults[[1]]$config <- config
    inferenceResults[[2]] <- NA
   }
  univQvals <- p.adjust(univPvals, method = "BH")
  univBY <- p.adjust(univPvals, method = "BY")
  nDiscoveries <- sum(univQvals < 0.05 & trueCoefs != 0)
  byDiscovers <- sum(univBY < 0.05 & trueCoefs != 0)
  bhFDR <- sum(univQvals < 0.05 & trueCoefs == 0) / max(sum(univQvals < 0.05), 1)
  byFDR <- sum(univBY < 0.05 & trueCoefs == 0) / max(sum(univBY < 0.05), 1)
  nNonNull <- sum(trueCoefs != 0)
  bh <- c(fdr = bhFDR, discoveries = nDiscoveries, nNonNull = nNonNull)
  by <- c(fdr = byFDR, discoveries = byDiscovers, nNonNull = nNonNull)
  
  inferenceResults[[length(inferenceResults)]] <- rbind(bh, by)
  return(inferenceResults)
}

# Running simulation --------
configA <- expand.grid(geneSizes = list(c(25, 55, 100)),
                       R = 0.02,
                       signal = c("weak"), 
                       sparsity = c(0.05, 0.1, 0.2),
                       rho = 0.8,
                       nSubjects = 10^4,
                       MAFthreshold = 0.1,
                       pvalThreshold = 0.05,
                       genomSize = 100,
                       pNullGenes = c(1 - 0.005),
                       seed = 151:200)
configB <- expand.grid(geneSizes = list(c(25, 55, 100)),
                       R = 0.02,
                       signal = c("weak"), 
                       sparsity = c(0.05, 0.1, 0.2),
                       rho = 0.8,
                       nSubjects = 10^4,
                       MAFthreshold = 0.1,
                       pvalThreshold = 0.05,
                       genomSize = 100,
                       pNullGenes = c(1 - 0.01),
                       seed = 151:200)
configurations <- rbind(configB, configA)
set.seed(3)
configurations <- configurations[order(runif(nrow(configurations))), ]

seed <- configurations[setting, ]$seed
genes <- configurations[setting, ]$genomSize
system.time(result <- runSim(configurations[setting, ], verbose = TRUE))
filename <- paste("results/skatGenome_F_", genes, 
                  "_genes_seed_", seed, 
                  "_setting_", setting, ".rds", sep = "")
saveRDS(result, file = filename)
