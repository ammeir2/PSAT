library(ggplot2)
library(dplyr)
library(PSAT)

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
      dat[, i] <- pnorm(rand[, i]) <= Maf[i]
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
    contrast <- sign(trueBeta)
    contrast[contrast == 0] <- sample(c(-1, 1), sum(contrast == 0), replace = TRUE)
    betacov <- solve(XtX)
    contsd <- as.numeric(sqrt(t(contrast) %*% betacov %*% contrast))
    threshold <- contsd * qnorm(c(pthreshold, 1- pthreshold))
    naiveBeta <- as.numeric(sampleLinearTest(1, trueBeta, betacov, contrast, threshold))
    Xy <- as.numeric(XtX %*% naiveBeta)

    # Inference -------------------
    sigma <- solve(XtX)
    time <- system.time(condFit <- mvnLinear(y = naiveBeta, sigma = sigma,
                                             contrast = contrast, threshold = threshold,
                                             selection = "two_sided",
                                             estimate_type = c("mle", "naive"),
                                             pvalue_type = "polyhedral",
                                             ci_type = c("polyhedral", "naive", "switch")))
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
                            t2type = t2type)
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
                              reps = 200,
                              betaType = 1,
                              noiseType = 0,
                              pthreshold = c(10^-3),
                              MAFthreshold = c(0.1),
                              t2type = c(2),
                              randRho = 0.8)
configurations <- configurations[sample(nrow(configurations)), ]
system.time(simResults <- apply(configurations, 1, run.simulation))
# save(simResults, file = "simulations/results/paper CI linear oct2G.Robj")

load(file = "simulations/results/paper CI linear oct2F.Robj")
simResultsF <- simResults
load(file = "simulations/results/paper CI linear oct2E.Robj")
simResultsE <- simResults
load(file = "simulations/results/paper CI linear oct2D.Robj")
simResultsD <- simResults
load(file = "simulations/results/paper CI linear oct2G.Robj")
simResultsG <- simResults

simResults <- c(simResultsD, simResultsE, simResultsG, simResultsF)

coverList <- list()
detList <- list()
sizeList <- list()
for(i in 1:length(simResults)) {
  iterResult <- simResults[[i]]
  cover <- t(sapply(iterResult, function(x) x$cover))
  det <- t(sapply(iterResult, function(x) x$det))
  cover <- data.frame(p = iterResult[[1]]$params[3],
                      signal = iterResult[[1]]$params[[5]],
                      t2 = iterResult[[1]]$params[[12]],
                      sparisty = iterResult[[1]]$params[[4]],
                      cover = cover, measure = "cover")
  det <- data.frame(p = iterResult[[1]]$params[3],
                    signal = iterResult[[1]]$params[[5]],
                    t2 = iterResult[[1]]$params[[12]],
                    sparisty = iterResult[[1]]$params[[4]],
                    det = det, measure = "signDet")

  names(cover) <- c("p", "signal", "t2", "sparsity", "naive", "switch", "poly", "measure")
  names(det) <- c("p", "signal", "t2", "sparsity", "naive", "switch", "poly", "measure")
  coverList[[i]] <- cover
  detList[[i]] <- det
}

library(dplyr)
library(reshape2)
library(ggplot2)
results <- rbind(do.call("rbind", coverList), do.call("rbind", detList))
results <- melt(results, id = c("p", "signal", "t2", "sparsity", "measure"))
names(results)[6:7] <- c("method", "error")
switch <- subset(results, method == "switch")
switch$method <- NULL
switch$method[switch$t2 == 1] <- "switch-sqrd"
switch$method[switch$t2 == 2] <- "switch-half"
switch$t2 <- NULL
results <- subset(results, method != "switch")
results$t2 <- NULL
results <- rbind(results, switch)
results <- summarize(group_by(results, p, signal, sparsity, method, measure),
                     sd = sd(error) / sqrt(length(error)),
                     error = mean(error))
results$error[results$signal == 0 & results$measure == "signDet"] <- 0
results$sd[results$signal == 0 & results$measure == "signDet"] <- 0
results$lci <- results$error - results$sd * 2
results$uci <- results$error + results$sd * 2

results$measureFacet <- "% Coverage Rate"
results$measureFacet[results$measure == "signDet"] <- "% Sign Determination"
results$sparseFacet <- paste("Sparsity:", results$sparsity)
results$intercept <- NA
results$intercept[results$measure == "cover"] <- 0.95
results$method <- as.character(results$method)
results$method[results$method == "switch-half"] <- "switch"
results$method[results$method == "poly"] <- "polyhedral"
results$method <- factor(results$method, levels = c("naive", "polyhedral", "switch"))

# pdf("figures/ciplotWpowerLinear_R.pdf",pagecentre=T, width=8, height=3.5 ,paper = "special")
results$R_squared <- results$signal / (1 + results$signal)
results$nonzeroCoef <- paste("Non-zero Coefs: ", results$sparsity)
ggplot(subset(results, sparsity != 9), aes(x = log2(R_squared), y = error, col = method, linetype = method)) +
  geom_line() + geom_point(aes(shape = method)) +
  facet_grid(measureFacet ~ nonzeroCoef, scales = "free_y") +
  geom_segment(aes(x = log2(R_squared), xend = log2(signal), y = lci, yend = uci)) +
  theme_bw() + geom_hline(aes(yintercept = intercept)) +
  #geom_hline(aes(yintercept = 1 - (1 - intercept) / 2), linetype = 2) +
  ylab("") + xlab("log2(R_squared)")
# dev.off()

# Coverage of zero vs. non-zero ---------------
# resultList <- c()
# for(i in 1:length(simResults)) {
#   polyCoverZero <- numeric(length(simResults[[i]]))
#   polyCoverAlt <- numeric(length(simResults[[i]]))
#   polyCoverAll <- numeric(length(simResults[[i]]))
#   naiveCoverZero <- numeric(length(simResults[[i]]))
#   naiveCoverAlt <- numeric(length(simResults[[i]]))
#   naiveCoverAll <- numeric(length(simResults[[i]]))
#   switchCoverZero <- numeric(length(simResults[[i]]))
#   switchCoverAlt <- numeric(length(simResults[[i]]))
#   switchCoverAll <- numeric(length(simResults[[i]]))
#   for(j in 1:length(simResults[[i]])) {
#     sparsity <- iterResult$params[4]
#     p <- iterResult$params[3]
#     iterResult <- simResults[[i]][[j]]
#     true <- iterResult$true
#     poly <- iterResult$tnCI
#     poly <- poly[, 1] < true & poly[, 2] > true
#     polyCoverZero[j] <- mean(poly[(sparsity + 1):p])
#     polyCoverAlt[j] <- mean(poly[1:sparsity])
#     polyCoverAll[j] <- mean(poly)
#
#     naive <- iterResult$naiveCI
#     naive <- naive[, 1] < true & naive[, 2] > true
#     naiveCoverZero[j] <- mean(naive[(sparsity + 1):p])
#     naiveCoverAlt[j] <- mean(naive[1:sparsity])
#     naiveCoverAll[j] <- mean(naive)
#
#     switch <- iterResult$condCI
#     switch <- switch[, 1] < true & switch[, 2] > true
#     switchCoverZero[j] <- mean(switch[(sparsity + 1):p])
#     switchCoverAlt[j] <- mean(switch[1:sparsity])
#     switchCoverAll[j] <- mean(switch)
#   }
#   all <- data.frame(n = iterResult$params[2],
#                     sparsity = iterResult$params[4],
#                     signal = iterResult$params[5],
#                     t2type = iterResult$params[12],
#                     poly = polyCoverAll,
#                     naive = naiveCoverAll,
#                     switch = switchCoverAll,
#                     subset = "all")
#   nonzero <- data.frame(n = iterResult$params[2],
#                         sparsity = iterResult$params[4],
#                         signal = iterResult$params[5],
#                         t2type = iterResult$params[12],
#                         poly = polyCoverAlt,
#                         naive = naiveCoverAlt,
#                         switch = switchCoverAlt,
#                         subset = "nonzero")
#   zero <- data.frame(n = iterResult$params[2],
#                      sparsity = iterResult$params[4],
#                      signal = iterResult$params[5],
#                      t2type = iterResult$params[12],
#                      poly = polyCoverZero,
#                      naive = naiveCoverZero,
#                      switch = switchCoverZero,
#                      subset = "zero")
#   iter <- rbind(all, zero, nonzero)
#   resultList[[i]] <- iter
# }
# results <- do.call("rbind", resultList)
# results <- melt(results, id = c("n", "sparsity", "signal", "t2type", "subset"))
# names(results)[6:7] <- c("method", "cover")
#
# sqrd <- subset(results, t2type == 2 & method == "switch")
# sqrd$method <- "switch-sqrd"
# half <- subset(results, t2type == 1 & method == "switch")
# half$method <- "switch-half"
# results <- subset(results, method != "switch")
# results <- rbind(results, half, sqrd)
# results$method <- factor(results$method, levels = c("naive", "poly",
#                                                     "switch-half", "switch-sqrd"))
#
#
# results <- summarize(group_by(results, n, sparsity, signal, subset, method),
#                      sd = sd(cover) / sqrt(length(cover)),
#                      cover = mean(cover))
#
# results$sparseFacet <- paste("Sparsity:", results$sparsity)
# results$subsetFacet <- paste("Coefficients:", results$subset)
#
# # pdf("figures/ciplotBySubset.pdf",pagecentre=T, width=8,height=3.5 ,paper = "special")
# ggplot(results, aes(x = signal, y = cover, col = method, linetype = method)) +
#   geom_line() +
#   facet_grid(sparseFacet ~ subsetFacet) +
#   theme_bw() +
#   geom_segment(aes(x = signal, xend = signal,
#                    y = cover + 2 * sd, yend  = cover - 2 * sd)) +
#   geom_hline(yintercept = 0.95) +
#   ylab("CI Coverage Rate") +
#   xlab("Siganl to Noise Ratio")
# # dev.off()


# Size of zero vs. non-zero ---------------
# resultList <- c()
# for(i in 1:length(simResults)) {
#   polyZero <- numeric(length(simResults[[i]]))
#   polyAlt <- numeric(length(simResults[[i]]))
#   polyAll <- numeric(length(simResults[[i]]))
#   naiveZero <- numeric(length(simResults[[i]]))
#   naiveAlt <- numeric(length(simResults[[i]]))
#   naiveAll <- numeric(length(simResults[[i]]))
#   switchZero <- numeric(length(simResults[[i]]))
#   switchAlt <- numeric(length(simResults[[i]]))
#   switchAll <- numeric(length(simResults[[i]]))
#   for(j in 1:length(simResults[[i]])) {
#     sparsity <- iterResult$params[4]
#     p <- iterResult$params[3]
#     iterResult <- simResults[[i]][[j]]
#     true <- iterResult$true
#     poly <- iterResult$tnCI
#     zeroind <- (sparsity + 1):p
#     altind <- 1:sparsity
#     polyZero[j] <- mean(poly[zeroind, 2] - poly[zeroind, 1])
#     polyAlt[j] <- mean(poly[altind, 2] - poly[altind, 1])
#     polyAll[j] <- mean(poly[, 2] - poly[, 1])
#
#     naive <- iterResult$naiveCI
#     naiveZero[j] <- mean(naive[zeroind, 2] - naive[zeroind, 1])
#     naiveAlt[j] <- mean(naive[altind, 2] - naive[altind, 1])
#     naiveAll[j] <- mean(naive[, 2] - naive[, 1])
#
#     switch <- iterResult$condCI
#     switchZero[j] <- mean(switch[zeroind, 2] - switch[zeroind, 1])
#     switchAlt[j] <- mean(switch[altind, 2] - switch[altind, 1])
#     switchAll[j] <- mean(switch[, 2] - switch[, 1])
#   }
#   all <- data.frame(n = iterResult$params[2],
#                     sparsity = iterResult$params[4],
#                     signal = iterResult$params[5],
#                     t2type = iterResult$params[12],
#                     poly = polyAll,
#                     naive = naiveAll,
#                     switch = switchAll,
#                     subset = "all")
#   nonzero <- data.frame(n = iterResult$params[2],
#                         sparsity = iterResult$params[4],
#                         signal = iterResult$params[5],
#                         t2type = iterResult$params[12],
#                         poly = polyAlt,
#                         naive = naiveAlt,
#                         switch = switchAlt,
#                         subset = "nonzero")
#   zero <- data.frame(n = iterResult$params[2],
#                      sparsity = iterResult$params[4],
#                      signal = iterResult$params[5],
#                      t2type = iterResult$params[12],
#                      poly = polyZero,
#                      naive = naiveZero,
#                      switch = switchZero,
#                      subset = "zero")
#   iter <- rbind(all, zero, nonzero)
#   resultList[[i]] <- iter
# }
# results <- do.call("rbind", resultList)
# results <- melt(results, id = c("n", "sparsity", "signal", "t2type", "subset"))
# names(results)[6:7] <- c("method", "size")
#
# sqrd <- subset(results, t2type == 2 & method == "switch")
# sqrd$method <- "switch-sqrd"
# half <- subset(results, t2type == 1 & method == "switch")
# half$method <- "switch-half"
# results <- subset(results, method != "switch")
# results <- rbind(results, half, sqrd)
# results$method <- factor(results$method, levels = c("naive", "poly",
#                                                     "switch-half", "switch-sqrd"))
#
#
# results <- summarize(group_by(results, n, sparsity, signal, subset, method),
#                      sd = sd(size) / sqrt(length(size)),
#                      size = mean(size))
#
# results$sparseFacet <- paste("Sparsity:", results$sparsity)
# results$subsetFacet <- paste("Coefficients:", results$subset)
#
# # pdf("figures/ciplotBySubset.pdf",pagecentre=T, width=8,height=3.5 ,paper = "special")
# ggplot(results, aes(x = signal, y = size, col = method, linetype = method)) +
#   geom_line(aes(shape = method)) +
#   facet_grid(sparseFacet ~ subsetFacet) +
#   theme_bw() +
#   geom_segment(aes(x = signal, xend = signal,
#                    y = size + 2 * sd, yend  = size - 2 * sd)) +
#   ylab("CI Size") +
#   xlab("Siganl to Noise Ratio")
# dev.off()
