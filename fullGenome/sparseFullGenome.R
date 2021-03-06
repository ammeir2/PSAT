# Information from cluster --- 
args <- commandArgs(TRUE)
eval(parse(text=args[[1]]))
setting <- as.numeric(setting)

# Helper functions ----
library(magrittr)
library(PSAT)
library(MASS)
library(Matrix)

checkInclusion <- function(ci, true) {
  inside <- 0
  for(i in 1:nrow(ci)) {
    if(ci[i, 1] < true[i] & ci[i, 2] > true[i]) {
      inside <- inside + 1
    }
  }
  
  return(inside / nrow(ci))
}

# Data generating functions ------
genNormVariants <- function(nSubjects, nVariants, rho, init = NULL) {
  if(is.null(init)) {
    init <- rnorm(nSubjects)
  }
  
  gene <- matrix(rnorm(nSubjects * nVariants), ncol = nVariants)
  for(j in 1:ncol(gene)) {
    if(j == 1) {
      gene[, 1] <- gene[, 1] + rho * init
    } else {
      gene[, j] <- gene[, j] + rho * gene[, j - 1]
    }
  }
  
  return(scale(gene))
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

genGene <- function(minGeneSize, maxGeneSize, nSubjects, rho, prevNorm,
                    minMAF, MAFthreshold, pNullGenes, 
                    sparseLevels, sparseProbs,
                    snrLevels, snrProbs,
                    gammaShape = 1, gammaRate = 50, seed = NULL) {
  if(!is.null(seed)) {
    set.seed(seed)
  }
  nVariants <- sample(minGeneSize:maxGeneSize, 1)
  normGene <- genNormVariants(2 * nSubjects, nVariants, rho, init = prevNorm)
  prevNorm <- normGene[, nVariants]
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
    mu <- as.numeric(X %*% trueCoef)
    mu <- mu / sd(mu) * snr
  } else {
    mu <- rep(0, nSubjects)
    nonzero <- 0
  }
  
  return(list(X = X, sparseX = sparseX, mu = mu, nonzero = nonzero))
}

# Simulation function --------
runSim <- function(config) {
  attach(config, warn.conflicts = FALSE)
  set.seed(seed)
  seeds <- sample.int(10^6, genomSize)

  # Sparsity ----
  sparseLevels <- c(1, 2, 4, 8)
  if(sparsity == "sparse") {
    sparseProbs <- c(4:1) / 10
  } else if(sparsity == "dense") {
    sparseProbs <- 1:4 / 10
  }
  
  # SNR levels -----
  snrLevels <- c(0.015, 0.03, 0.45, 0.6) 
  if(signal == "weak") {
    snrProbs <- 4:1 / 10
  } else if(signal == "strong") {
    snrProbs <- 1:4 / 10
  }
  
  # computing expected value ------
  yMean <- rep(0, nSubjects)
  
  print("Generating response")
  prevNorm <- NULL
  sparseMats <- list(genomSize, mode = "list")
  pb <- txtProgressBar(min = 0, max = genomSize, style = 3)
  for(g in 1:genomSize) {
    # Generating data
    dat <- genGene(minGeneSize, maxGeneSize, nSubjects, rho, prevNorm,
                   minMAF, MAFthreshold, pNullGenes, 
                   sparseLevels, sparseProbs,
                   snrLevels, snrProbs,
                   gammaShape = 1, gammaRate = 50,
                   seed = seeds[g])  
    sparseMats[[g]] <- dat$sparseX
    prevNorm <- dat$X[, ncol(dat$X)]
    yMean <- yMean + dat$mu
    setTxtProgressBar(pb, g)
  }
  rm(dat)
  close(pb)
  yMean <- yMean / sd(yMean) * totSNR
  y <- rnorm(nSubjects, mean = yMean, sd = 1)
  yMean <- yMean - mean(y)
  y <- y - mean(y)
  yvar <- var(y)
  
  print("Aggregate Testing")
  slot <- 1
  set.seed(seed)
  prevNorm <- NULL
  selectedGenes <- vector(genomSize, mode = "list")
  pb <- txtProgressBar(min = 0, max = genomSize, style = 3)
  for(g in 1:genomSize) {
    # Getting sparse Matrix 
    X <- sparseMats[[g]] %>% as.matrix() %>% scale()
    sparseMats[[g]] <- 0 # saving memory

    # Aggregate testing -------
    XtX <- t(X) %*% X
    XtXinv <- ginv(XtX)
    suffStat <- t(X) %*% y
    naiveCoef <- XtXinv %*% suffStat
    coefCov <- yvar * XtXinv
    waldMat <- XtX / as.numeric(var(y - X %*% naiveCoef))
    waldStat <- as.numeric(t(naiveCoef) %*% waldMat %*% naiveCoef)
    critVal <- qchisq(pvalThreshold / genomSize, df = ncol(X), lower.tail = FALSE)
    trueCoef <- XtXinv %*% t(X) %*% yMean
    naiveSD <- sqrt(diag(coefCov)) * sd(y - as.numeric(X %*% naiveCoef)) / sd(y)
    # print(c(g, ncol(X), waldStat, critVal, sum(abs(trueCoef)), dat$nonzero))
    if(waldStat > critVal) {
      selectedGenes[[slot]] <- list(obs = naiveCoef, 
                                    cov = coefCov, 
                                    naiveSD = naiveSD,
                                    true = trueCoef,
                                    waldMat = waldMat,
                                    threshold = critVal)
      slot <- slot + 1
    }
    setTxtProgressBar(pb, g)
  }
  rm(X)
  close(pb)
  
  if(slot == 1) {
    result <- list()
    result[[1]] <- list()
    result[[1]]$config <- config
    print("No significant genes found!")
    return(result)
  }
  
  
  # Conducting Inference
  selectedGenes <- selectedGenes[1:(slot - 1)]
  inferenceResults <- vector(length(selectedGenes), mode = "list")
  totVariants <- sapply(selectedGenes, function(x) ncol(x$cov)) %>% sum()
  print("Conducting Inference")
  pb <- txtProgressBar(min = 0, max = length(selectedGenes), style = 3)
  for(g in 1:length(selectedGenes)) {
    gene <- selectedGenes[[g]]
    naive <- gene$obs %>% as.numeric()
    true <- gene$true %>% as.numeric()
    sds <- sqrt(diag(gene$cov))
    naiveSD <- gene$naiveSD
    truePvals <- 2 * pnorm(-abs(true /sds))
    empNonzero <- sum(truePvals < 0.05)
    
    psatFit <- mvnQuadratic(naive, gene$cov, testMat = gene$waldMat, 
                            threshold = gene$threshold, 
                            contrasts = diag(length(naive)),
                            estimate_type = "naive",
                            pvalue_type = c("hybrid", "polyhedral"), 
                            ci_type = c("switch", "polyhedral"),
                            verbose = FALSE)
    # Cover Rate
    switchCover <- checkInclusion(psatFit$switchCI, true)
    polyCover <- checkInclusion(psatFit$polyCI, true)
    bbQuantile <- length(selectedGenes) * 0.05 / 2 / genomSize
    bbCritVal <- qnorm(1 - bbQuantile)
    bbCI <- naive + cbind(-naiveSD * bbCritVal, naiveSD * bbCritVal)
    bbCover <- checkInclusion(bbCI, true)
    
    # Avg Power 
    hybridAvgPower <- min(sum(p.adjust(psatFit$hybridPval, method = "BH") < 0.05) / empNonzero, 1)
    polyAvgPower <- min(sum(p.adjust(psatFit$polyPval, method = "BH") < 0.05) / empNonzero, 1)
    bbPval <- 2 * pnorm(-abs(naive / naiveSD))
    bbAvgPower <- min(sum(p.adjust(bbPval, method = "BH") < bbQuantile) / empNonzero, 1)
    
    # Reporting
    geneResult <- list()
    geneResult$config <- config
    geneResult$cover <- c(switch = switchCover, poly = polyCover, bb = bbCover)
    geneResult$power <- c(hybrid = hybridAvgPower, poly = polyAvgPower, bbAvgPower = bbAvgPower)
    inferenceResults[[g]] <- geneResult

    setTxtProgressBar(pb, g)
  }
  close(pb)
  
  return(inferenceResults)
}

# Running simulation --------
configurations <- expand.grid(geneSizes = c(10, 50, 100),
                              totSNR = c(2, 1, 0.5, 0.25),
                              signal = c("weak"), 
                              sparsity = c("sparse", "dense"),
                              rho = 0.8,
                              nSubjects = 10^4,
                              MAFthreshold = 0.1,
                              pvalThreshold = 0.05,
                              genomSize = 20000,
                              pNullGenes = c(1 - 0.025, 1 - 0.0025),
                              seed = 101:200)
settings <- which(configurations$totSNR %in% c(2, 0.25))
result <- runSim(configurations[setting, ])
filename <- paste("results/fullGenome_D_20000genes_setting", setting, "B.rds", sep = "")
saveRDS(result, file = filename)
