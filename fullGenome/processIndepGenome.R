library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)

compAvgFDR <- function(x, levels = c(0.01, 0.05, 0.1, 0.15, 0.2)) {
  tab <- x$FDR
  bbFraction <- tab$bbFraction[1]
  result <- matrix(ncol = 5, nrow = length(levels))
  for(i in 1:length(levels)) {
    level <- levels[i]
    qvals <- apply(tab[, 1:3], 2, p.adjust, method = "BH")
    hybrid <- sum(qvals[, 1] < level & tab$true == 0) / sum(qvals[, 1] < level)
    poly <- sum(qvals[, 2] < level & tab$true == 0) / sum(qvals[, 2] < level)
    bb <- sum(qvals[, 3] < level * bbFraction & tab$true == 0) / sum(qvals[, 3] < level * bbFraction)
    naive <- sum(qvals[, 3] < level & tab$true == 0) / sum(qvals[, 3] < level)
    result[i, ] <- c(hybrid = hybrid, poly = poly, bb = bb, naive = naive, level = level)
  }
  result[is.nan(result)] <- 0
  result <- data.frame(result)
  names(result) <- c("hybrid", "poly", "bb", "naive", "level")
  result <- suppressWarnings(cbind(data.frame(x$config), result))
  return(result)
}

compOverallPower <- function(x, level = 0.05) {
  len <- length(x)
  res <- lapply(x[-len], function(y) y$FDR) %>% do.call("rbind", .)
  res[, 1:3] <- apply(res[, 1:3], 2, p.adjust, method = "BH")
  hybridFDR <- sum(res$hybrid < 0.05 & res$true == 0) / sum(res$hybrid < 0.05)
  polyFDR <- sum(res$poly < 0.05 & res$true == 0) / sum(res$poly < 0.05)
  naiveFDR <- sum(res$naive < 0.05 & res$true == 0) / sum(res$naive < 0.05)
  res <- subset(res, true != 0)
  hybrid <- sum(res$hybrid < 0.05)
  polyhedral <- sum(res$poly < 0.05)
  bh <- x[[length(x)]]
  bhFDR <- as.numeric(bh[1, 1])
  byFDR = as.numeric(bh[2, 1])
  power <- c(hybrid = hybrid, polyhedral = polyhedral, bh = bh[1, 2], by = bh[2, 2]) / bh[1, 3]
  fdr <- c(hybrid = hybridFDR, poly = polyFDR, naive = naiveFDR, bh = bh[1, 1], by = bh[2, 1])
  fdr[is.nan(fdr)] <- 0
  return(list(fdr = fdr, power = power))
}

filenames <- as.list(dir(path = 'fullGenome/results', pattern="bhGenome_D_*"))
filenames <- lapply(filenames, function(x) paste0('fullGenome/results/', x))
filenames <- as.character(filenames)
results <- vector(length(filenames), mode = "list")
pb <- txtProgressBar(min = 0, max = length(filenames), style = 3)
for(i in 1:length(filenames)) {
  result <- readRDS(filenames[[i]])
  if(length(result[[1]]) == 1) {
    next
  }
  len <- length(result)
  dat <- lapply(result[-len], function(x) x$config) %>% do.call("rbind", .)
  power <- lapply(result[-len], function(x) x$power) %>% do.call("rbind", .)
  selectedNonNull <- sapply(result[-len], function(x) sum(x$FDR$true != 0))
  nTrueReject <- apply(power, 2, function(x) x * selectedNonNull) %>% 
    matrix(ncol = 3) %>% colSums()
  totpower <- nTrueReject / result[[len]][1, 3]
  names(totpower) <- names(result[[1]]$power)
  cover <- lapply(result[-len], function(x) x$cover) %>% do.call("rbind", .)
  avgFDR <- lapply(result[-len], compAvgFDR) %>% do.call("rbind", .)
  overallPower <- compOverallPower(result)
  overallFDR <- overallPower$fdr
  overallPower <- overallPower$power
  # print(c(i, nrow(dat), nrow(power)))
  if(nrow(dat) != nrow(power)) {
    print("hello")
  }
  results[[i]] <- list(dat = dat, power = power, cover = cover, 
                       totpower = totpower,
                       fdr = avgFDR,
                       overallPower = overallPower,
                       overallFDR = overallFDR)
  setTxtProgressBar(pb, i)
}
close(pb)
results <- results[!sapply(results, is.null)]

# Overall power -------
config <- do.call("rbind", lapply(results, function(x) x$dat[1, ]))
bhCompare <- do.call("rbind", lapply(results, function(x) x$overallPower))
bhCompare <- cbind(config, bhCompare)
bhCompare <- melt(bhCompare, id = colnames(bhCompare)[-c(12:15)])
names(bhCompare)[12:13] <- c("method", "power")
bhCompare <- group_by(bhCompare, totSNR, sparsity, pNullGenes, method) %>% 
  summarize(powerSD = sd(power) / sqrt(length(power)), 
            power = mean(power)) %>% 
  data.frame()
# bhCompare <- subset(bhCompare, !(log2(totSNR) == -2 & pNullGenes == 0.975))
ggplot(bhCompare, aes(x = log2(totSNR), y = power, col = method, linetype = method, shape = method)) + 
  geom_point() + geom_line() +
  geom_segment(aes(y = power - 2 * powerSD, yend = power + 2 * powerSD, xend = log2(totSNR))) + 
  theme_bw() + 
  facet_grid(sparsity ~ pNullGenes, scales = "free_x", labeller = "label_both")
ggsave(file = "figures/5000genesBHpower.pdf", width = 10, height = 5)

# Overall FDR ----------
config <- do.call("rbind", lapply(results, function(x) x$dat[1, ]))
totFDR <- do.call("rbind", lapply(results, function(x) x$overallFDR))
totFDR[is.nan(totFDR)] <- 0
totFDR <- cbind(config, totFDR)
totFDR <- gather(totFDR, key = method, value = fdr, hybrid, poly, naive, bh, by)
totFDR <- group_by(totFDR, totSNR, sparsity, pNullGenes, method) %>% 
  summarize(fdrSD = sd(fdr) / sqrt(length(fdr)), fdr = mean(fdr))
# totFDR <- subset(totFDR, !(log2(totSNR) == -2 & pNullGenes == 0.975))
ggplot(totFDR, aes(x = log2(totSNR), y = fdr, col = method, linetype = method, shape = method)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid(sparsity ~ pNullGenes, labeller = "label_both", scales = "free_x") + 
  geom_segment(aes(xend = log2(totSNR), 
                   yend = pmin(fdr + 2 * fdrSD, 1), y = pmax(fdr - 2 * fdrSD, 0))) + 
  geom_hline(yintercept = 0.05) +
  ggtitle("Overall FDR")
ggsave(file = "figures/5000genesOverallFDR.pdf", width = 10, height = 5)

# Power ----------
power <- do.call("rbind", lapply(results, function(x) x$power))
config <- do.call("rbind", lapply(results, function(x) x$dat))
power <- cbind(config, power)
nGenes <- group_by(config, totSNR, signal, sparsity, pNullGenes) %>% 
  summarize(discoveries = length(seed) / length(unique(seed)))
power <- melt(power, id = colnames(power)[-(12:14)])
names(power)[12:13] <- c("method", "power")
power <- subset(power, !is.nan(power))
# power <- subset(power, !(log2(totSNR) == -2 & pNullGenes == 0.975))
power <- group_by(power, totSNR, sparsity, pNullGenes, method, seed) %>% 
  summarize(power = mean(power)) %>% 
  group_by(totSNR, sparsity, pNullGenes, method) %>% 
  summarize(powerSD = sd(power, na.rm = TRUE) / sqrt(sum(!is.na(power))), 
            power = mean(power, na.rm = TRUE))

ggplot(power, aes(x = log2(totSNR), y = power, col = method, linetype = method, shape = method)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid(sparsity ~ pNullGenes, labeller = "label_both", scales = "free_x") + 
  geom_segment(aes(xend = log2(totSNR), 
                   yend = pmin(power + 2 * powerSD, 1), y = pmax(power - 2 * powerSD, 0))) + 
  ggtitle("Power of AVG FDR (Over Selected)")
ggsave(file = "figures/5000genesPower_overSelected.pdf", width = 10, height = 5)

# FDR ----------
fdr <- do.call("rbind", lapply(results, function(x) x$fdr))
fdr <- melt(fdr, id = colnames(fdr)[-(12:15)])
names(fdr)[13:14] <- c("method", "fdr")
fdr <- subset(fdr, !is.nan(fdr))
fdr <- group_by(fdr, totSNR, sparsity, pNullGenes, method, seed, level) %>% 
  summarize(fdr = mean(fdr)) %>% 
  group_by(totSNR, sparsity, pNullGenes, method, level) %>% 
  summarize(fdrSD = sd(fdr, na.rm = TRUE) / sqrt(sum(!is.na(fdr))), 
            fdr = mean(fdr, na.rm = TRUE))
fdr <- fdr[order(fdr$totSNR), ]
# fdr <- subset(fdr, !(log2(totSNR) == -2 & pNullGenes == 0.975))
fdrDense <- subset(fdr, pNullGenes == 0.975)
fdrDense$snr_level <- plyr::mapvalues(fdrDense$totSNR, from = sort(unique(fdrDense$totSNR)), to = 1:length(unique(fdrDense$totSNR)))
fdrSparse <- subset(fdr, pNullGenes == .9975)
fdrSparse$snr_level <- plyr::mapvalues(fdrSparse$totSNR, from = sort(unique(fdrSparse$totSNR)), to = 1:length(unique(fdrSparse$totSNR)))
fdr <- rbind(fdrDense, fdrSparse)

normQuant <- qnorm(1 - 0.05 / nrow(fdr))
ggplot(fdr, aes(x = level, y = fdr, col = method, linetype = method, shape = method)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid(snr_level ~ pNullGenes + sparsity, labeller = "label_both") + 
  geom_segment(aes(xend = level, 
                   yend = pmin(fdr + normQuant * fdrSD, 1), 
                   y = pmax(fdr - normQuant * fdrSD, 0))) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  ggtitle("Average FDR")
ggsave(file = "figures/5000genesfdr.pdf", width = 7, height = 6)

# Power ----------
power <- do.call("rbind", lapply(results, function(x) x$totpower))
config <- do.call("rbind", lapply(results, function(x) x$dat[1, ]))
power <- cbind(config, power)
nGenes <- group_by(config, totSNR, signal, sparsity, pNullGenes) %>% 
  summarize(discoveries = length(seed) / length(unique(seed)))
power <- melt(power, id = colnames(power)[-(12:14)])
names(power)[12:13] <- c("method", "power")
power <- subset(power, !is.nan(power))
# power <- subset(power, !(log2(totSNR) == -2 & pNullGenes == 0.975))
power <- group_by(power, totSNR, sparsity, pNullGenes, method, seed) %>% 
  summarize(power = mean(power)) %>% 
  group_by(totSNR, sparsity, pNullGenes, method) %>% 
  summarize(powerSD = sd(power, na.rm = TRUE) / sqrt(sum(!is.na(power))), 
            power = mean(power, na.rm = TRUE))

ggplot(power, aes(x = log2(totSNR), y = power, col = method, linetype = method, shape = method)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid(sparsity ~ pNullGenes, labeller = "label_both", scales = "free_x") + 
  geom_segment(aes(xend = log2(totSNR), 
                   yend = pmin(power + 2 * powerSD, 1), y = pmax(power - 2 * powerSD, 0))) + 
  ggtitle("Power of AVG FDR (Overall)")
ggsave(file = "figures/5000genesPower_overall.pdf", width = 10, height = 5)
