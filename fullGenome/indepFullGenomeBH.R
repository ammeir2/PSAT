# Information from cluster ---
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
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
                    minMAF, MAFthreshold, pNullGenes, 
                    sparseLevels, sparseProbs,
                    snrLevels, snrProbs,
                    gammaShape = 1, gammaRate = 50, seed = NULL) {
  if(!is.null(seed)) {
    set.seed(seed)
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
  nonNull <- runif(1) <= 1 - pNullGenes
  trueCoef <- rep(0, nVariants)
  if(nonNull) {
    nonzero <- sample(sparseLevels, 1, FALSE, sparseProbs)
    snr <- sample(snrLevels, 1, FALSE, snrProbs)
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
  set.seed(seed)
  seeds <- sample.int(10^6, genomSize)
  sqrtMats <- lapply(unlist(geneSizes), function(size) {
    covmat <- rho^as.matrix(dist(1:size))
    sqrtmat <- sqrtm(covmat)
    return(sqrtmat)
  })

  # Sparsity ----
  sparseLevels <- c(1, 2, 4, 8)
  if(sparsity == "sparse") {
    sparseProbs <- c(4:1) / 10
  } else if(sparsity == "dense") {
    sparseProbs <- 1:4 / 10
  } else if(sparsity == "densest") {
    sparseProbs <- c(0, 0, 0, 1)
  }
  
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
  if(verbose) pb <- txtProgressBar(min = 0, max = genomSize, style = 3)
  totVariants <- 0  
  for(g in 1:genomSize) {
    # Generating data
    dat <- genGene(sqrtMats, nSubjects,
                   minMAF, MAFthreshold, pNullGenes, 
                   sparseLevels, sparseProbs,
                   snrLevels, snrProbs,
                   gammaShape = 1, gammaRate = 50,
                   seed = seeds[g])  
    sparseMats[[g]] <- dat$sparseX
    yMean <- yMean + dat$mu
    totVariants <- totVariants + ncol(dat$X)
    if(verbose) setTxtProgressBar(pb, g)
  }
  rm(dat)
  if(verbose) close(pb)
  yMean <- yMean / sd(yMean) * totSNR
  y <- rnorm(nSubjects, mean = yMean, sd = 1)
  yMean <- yMean - mean(y)
  y <- y - mean(y)
  yvar <- var(y)
  
  print("Aggregate Testing")
  slot <- 1
  set.seed(seed)
  selectedGenes <- vector(genomSize, mode = "list")
  if(verbose) pb <- txtProgressBar(min = 0, max = genomSize, style = 3)
  univPvals <- numeric(totVariants)
  trueCoefs <- numeric(totVariants)
  univInd <- 1
  for(g in 1:genomSize) {
    dat <- genGene(sqrtMats, nSubjects,
                   minMAF, MAFthreshold, pNullGenes, 
                   sparseLevels, sparseProbs,
                   snrLevels, snrProbs,
                   gammaShape = 1, gammaRate = 50,
                   seed = seeds[g])  
    X <- dat$X
    trueCoef <- dat$coef
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
                                    waldPval = waldPval)
      slot <- slot + 1
    }
    if(verbose) setTxtProgressBar(pb, g)
  }
  rm(X)
  if(verbose) close(pb)
  
  if(slot == 1) {
    result <- list()
    result[[1]] <- list()
    result[[1]]$config <- config
    print("No significant genes found!")
    return(result)
  }
  
  # selection with BH
  selectedGenes <- selectedGenes[1:(slot - 1)]
  aggPvals <- sapply(selectedGenes, function(x) x$waldPval)
  aggQvals <- p.adjust(aggPvals, method = "BH", n = genomSize)
  keep <- aggQvals < pvalThreshold
  bhFrac <- sum(keep) / genomSize
  selectedGenes <- selectedGenes[keep]
  for(i in 1:length(selectedGenes)) {
    geneDF <- length(selectedGenes[[i]]$obs)
    selectedGenes[[i]]$threshold <- qchisq(1 - pvalThreshold * bhFrac, df = geneDF)
  }
  
  # Conducting Inference
  inferenceResults <- vector(length(selectedGenes) + 1, mode = "list")
  totVariants <- sapply(selectedGenes, function(x) ncol(x$cov)) %>% sum()
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
    
    control <- psatControl(nSamples = 10^4)
    psatFit <- mvnQuadratic(naive, gene$cov, testMat = gene$waldMat, 
                            threshold = gene$threshold, 
                            contrasts = diag(length(naive)),
                            estimate_type = "naive",
                            pvalue_type = c("hybrid", "polyhedral"), 
                            ci_type = c("switch", "polyhedral"),
                            verbose = FALSE,
                            control = control)
    # Cover Rate
    switchCover <- checkInclusion(psatFit$switchCI, trueProj)
    polyCover <- checkInclusion(psatFit$polyCI, trueProj)
    bbQuantile <- length(selectedGenes) * 0.05 / 2 / genomSize
    bbCritVal <- qnorm(1 - bbQuantile)
    bbCI <- naive + cbind(-naiveSD * bbCritVal, naiveSD * bbCritVal)
    bbCover <- checkInclusion(bbCI, trueProj)
    
    # FDR
    fdrFrame <- data.frame(hybrid = getPval(psatFit, type = "hybrid"), 
                           poly = getPval(psatFit, type = "polyhedral"), 
                           naive = getPval(psatFit, type = "naive"), 
                           bbFraction = length(selectedGenes) / genomSize, 
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
    geneResult$config <- config
    geneResult$cover <- c(switch = switchCover, poly = polyCover, bb = bbCover)
    geneResult$power <- c(hybrid = hybridPower, poly = polyPower, bbAvgPower = bbPower)
    geneResult$FDR <- fdrFrame
    inferenceResults[[g]] <- geneResult

    if(verbose) setTxtProgressBar(pb, g)
  }
  if(verbose) close(pb)
  
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
configA <- expand.grid(geneSizes = list(c(10, 55, 100)),
                       totSNR = 2^seq(from = -3, to = -1),
                       # totSNR = 1,
                       signal = c("weak"), 
                       sparsity = c("sparse", "dense", "densest"),
                       rho = 0.8,
                       nSubjects = 10^4,
                       MAFthreshold = 0.1,
                       pvalThreshold = 0.05,
                       genomSize = 5000,
                       pNullGenes = c(1 - 0.0025),
                       seed = 81:120)
configB <- expand.grid(geneSizes = list(c(10, 55, 100)),
                       totSNR = 2^seq(from = -1, to = 1),
                       signal = c("weak"), 
                       sparsity = c("sparse", "dense", "densest"),
                       rho = 0.8,
                       nSubjects = 10^4,
                       MAFthreshold = 0.1,
                       pvalThreshold = 0.05,
                       genomSize = 5000,
                       pNullGenes = c(1 - 0.025),
                       seed = 81:120)
configurations <- rbind(configB, configA)
# set.seed(1)
# configurations <- configurations[order(runif(nrow(configurations))), ]

seed <- configurations[setting, ]$seed
genes <- configurations[setting, ]$genomSize
system.time(result <- runSim(configurations[setting, ], verbose = TRUE))
filename <- paste("results/bhGenome_D_", genes, 
                  "genes_seed", seed, 
                  "_setting", setting, ".rds", sep = "")
saveRDS(result, file = filename)
