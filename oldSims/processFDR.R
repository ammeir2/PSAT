library(magrittr)
library(ggplot2)
library(dplyr)
library(reshape2)

filenames <- as.list(dir(path = 'oldSims/results', pattern="fdr_sim_A_*"))
filenames <- lapply(filenames, function(x) paste0('oldSims/results/', x))
filenames <- as.character(filenames)
results <- list(length(filenames))
for(i in 1:length(filenames)) {
  results[[i]] <- readRDS(filenames[[i]])
}

results <- do.call("c", do.call("c", results))
processResults <- function(res, levels = c(0.01, 0.05, 0.1, 0.15, 0.2)) {
  par <- res$params
  sparsity <- par[[4]]
  p <- par[[3]]
  parnames <- names(par)
  par <- data.frame(matrix(par, nrow = 1))
  names(par) <- parnames
  par <- data.frame(par, level = levels)

  poly <- p.adjust(c(res$polyAlt, res$polyNull), method = "BH")
  poly <- sapply(levels, function(level) sum(poly[(sparsity + 1):p] < level) / sum(poly < level))
  poly[is.nan(poly)] <- 0
  naive <- p.adjust(c(res$naiveAltPvals, res$naiveNullPvals), method = "BH")
  naive <- sapply(levels, function(level) sum(naive[(sparsity + 1):p] < level) / sum(naive < level))
  naive[is.nan(naive)] <- 0
  null <- p.adjust(c(res$conditionalAltPvals, res$conditionalNullPvals), method = "BH")
  null <- sapply(levels, function(level) sum(null[(sparsity + 1):p] < level) / sum(null < level))
  null[is.nan(null)] <- 0
  hybrid <- p.adjust(c(res$altNotScore, res$nullNotScore), method = "BH")
  hybrid <- sapply(levels, function(level) sum(hybrid[(sparsity + 1):p] < level) / sum(hybrid < level))
  hybrid[is.nan(hybrid)] <- 0
  par$polyhedral <- poly
  par$naive <- naive
  par$hybrid <- hybrid
  par$null <- null
  return(par)
}
processResultsBonf <- function(res, levels = c(0.01, 0.05, 0.1, 0.15, 0.2)) {
  par <- res$params
  sparsity <- par[[4]]
  p <- par[[3]]
  parnames <- names(par)
  par <- data.frame(matrix(par, nrow = 1))
  names(par) <- parnames
  par <- data.frame(par, level = levels)

  poly <- p.adjust(c(res$polyAlt, res$polyNull), method = "bonferroni")
  poly <- sapply(levels, function(level) any(poly[(sparsity + 1):p] < level))
  poly[is.nan(poly)] <- 0
  naive <- p.adjust(c(res$naiveAltPvals, res$naiveNullPvals), method = "bonferroni")
  naive <- sapply(levels, function(level) any(naive[(sparsity + 1):p] < level))
  naive[is.nan(naive)] <- 0
  null <- p.adjust(c(res$conditionalAltPvals, res$conditionalNullPvals), method = "bonferroni")
  null <- sapply(levels, function(level) any(null[(sparsity + 1):p] < level))
  null[is.nan(null)] <- 0
  hybrid <- p.adjust(c(res$altNotScore, res$nullNotScore), method = "bonferroni")
  hybrid <- sapply(levels, function(level) any(hybrid[(sparsity + 1):p] < level))
  hybrid[is.nan(hybrid)] <- 0
  par$polyhedral <- poly
  par$naive <- naive
  par$hybrid <- hybrid
  par$null <- null
  return(par)
}

results <- lapply(results, processResults)
results <- data.table::rbindlist(results)

library(dplyr)
library(ggplot2)
library(reshape2)
results$m <- NULL
results$selectProb <- NULL
results$chisqPval <- NULL
results$coefType <- NULL
results$MAFthreshold <- NULL
results$pthreshold <- NULL
results <- melt(results, id = c("n", "p", "sparsity", "signal", "level", "noiseType", "R_squared"))
names(results)[8:9] <- c("method", "fdr")
results <- summarize(group_by(results, n, p, sparsity, signal, method, level, noiseType),
                     R_squared = mean(R_squared),
                     sdfdr = sd(fdr) / sqrt(length(fdr)),
                     fdr = mean(fdr))
results$lci <- results$fdr - results$sdfdr * 2
results$uci <- results$fdr + results$sdfdr * 2
results$noisePlot <- "Gaussian"
results$noisePlot[results$noiseType == 2] <- "Laplace"
results$noisePlot[results$noiseType == 3] <- "Uniform"
results$noisePlot <- factor(results$noisePlot,
                            levels = c("Gaussian", "Laplace", "Uniform"))
# results$snr <- paste("snr:", results$signal)
results$R_squared <- paste("R_squared:", round(results$R_squared, 4))
results$method <- as.character(results$method)
results$method[results$method == "null"] <- "global-null"
results$method <- factor(results$method, levels = c("naive", "polyhedral", "global-null", "hybrid"))

# pdf("figures/fdrplot_B.pdf",pagecentre=T, width=8,height=3.5 ,paper = "special")
ggplot(subset(results, signal != 0.064),
       aes(x = level, y = fdr, col = method,
                    linetype = method, shape = method)) +
  geom_line() + geom_point() + theme_bw() +
  ylim(0, 1) +
  geom_segment(aes(xend = level, y = lci, yend = uci)) +
  facet_grid(R_squared ~ noisePlot) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) + xlab("Nominal FDR") +
  ylab("Empirical FDR")
# dev.off()




