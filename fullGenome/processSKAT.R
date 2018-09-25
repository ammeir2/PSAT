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
  res <- lapply(x[-len], function(y) cbind(y$FDR, test = y$config[["test"]])) %>%
    do.call("rbind", .)
  wald <- subset(res, test == "wald")
  skat <- subset(res, test == "skat")
  if(any(res$test == "skat")) {
    skat[, 1:3] <- apply(skat[, 1:3], 2, p.adjust, method = "BH")
  }
  if(any(res$test == "wald")) {
    wald[, 1:3] <- apply(wald[, 1:3], 2, p.adjust, method = "BH")
  }
  res <- rbind(skat, wald)
  res <- gather(res, key = "method", value = "qval", hybrid, poly, naive)
  res <- group_by(res, method, test) %>% 
    summarize(fdr = sum(qval < 0.05 & true == 0) / max(sum(qval < 0.05), 1),
              power = sum(qval < 0.05 & true != 0))
  bh <- x[[length(x)]]
  res <- data.frame(res)
  res <- rbind(res, 
               data.frame(method = c("bh", "by"), test = "all", fdr = bh[1:2, 1], power = bh[1:2, 2]))
  res$power <- res$power / bh[1, 3]
  return(res)
}

filenames <- as.list(dir(path = 'fullGenome/results', pattern="skatGenome_I_*"))
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
  if(len == 2) {
    # browser()
  }
  dat <- lapply(result[-len], function(x) x$config) %>% do.call("rbind", .)
  power <- lapply(result[-len], function(x) data.frame(matrix(x$power, nrow = 1), x$config$test)) %>% do.call("rbind", .)
  names(power) <- c("hybrid", "polyhedral", "bbAvgPower", "test")
  selectedNonNull <- sapply(result[-len], function(x) sum(x$FDR$true != 0))
  totpower <- power
  totpower[, 1:3] <- apply(totpower[, 1:3], 2, function(x) x * selectedNonNull) %>% 
    matrix(ncol = 3)
  totpower <- by(totpower[, 1:3], totpower$test, colSums) %>% as.list() %>% do.call("rbind", .)
  totpower <- data.frame(totpower / result[[len]][1, 3])
  totpower$test <- rownames(totpower)
  if(!("wald" %in% totpower$test)) {
    totpower <- rbind(totpower, data.frame(hybrid = 0, polyhedral = 0, bbAvgPower = 0, test = "wald"))
  }
  if(!("skat" %in% totpower$test)) {
    totpower <- rbind(totpower, data.frame(hybrid = 0, polyhedral = 0, bbAvgPower = 0, test = "skat"))
  }
  cover <- lapply(result[-len], function(x) x$cover) %>% do.call("rbind", .)
  avgFDR <- lapply(result[-len], compAvgFDR) %>% do.call("rbind", .)
  
  overall <- compOverallPower(result)
  if(!("wald" %in% unique(overall$test))) {
    waldres <- data.frame(method = c("hybrid", "naive", "poly"), 
                          test = "wald", fdr = 0, power = 0)
    overall <- rbind(overall, waldres)
  }
  if(!("skat" %in% overall$test)) {
    skatrest <- data.frame(method = c("hybrid", "naive", "poly"), 
                          test = "skat", fdr = 0, power = 0)
    overall <- rbind(overall, skatrest)
  }
  
  # print(c(i, nrow(dat), nrow(power)))
  if(nrow(dat) != nrow(power)) {
    print("hello")
  }
  results[[i]] <- list(dat = dat, power = power, cover = cover, 
                       totpower = totpower,
                       fdr = avgFDR,
                       overall = overall)
  setTxtProgressBar(pb, i)
}
close(pb)
results <- results[!sapply(results, is.null)]
# saveRDS(results, file = "fullGenome/results/genome_processed_I.rds")
results <- readRDS(file = "fullGenome/results/genome_processed_I.rds")

# Overall power -------
config <- lapply(results, function(x) dplyr::select(x$dat[1, ], -test, - geneSizes))
bhCompare <- lapply(results, function(x) x$overall)
res <- vector(length(bhCompare), mode = "list")
for(i in 1:length(res)) {
  res[[i]] <- cbind(config[[i]], bhCompare[[i]])
}
bhCompare <- do.call("rbind", res)
bhCompare <- group_by(bhCompare, R, sparsity, pNullGenes, method, test) %>% 
  summarize(powerSD = sd(power) / sqrt(length(power)), 
            power = mean(power),
            fdrSD = sd(fdr) / sqrt(length(fdr)), 
            fdr = mean(fdr)) %>% 
  data.frame()
# bhCompare <- subset(bhCompare, !(log2(totSNR) == -2 & pNullGenes == 0.975))
bhCompare$testMethod <- interaction(bhCompare$test, bhCompare$method)
bhCompare$percNonNull <- paste("Non-zero variants: ", bhCompare$sparsity * 100, "%", sep = "")
bhCompare$percGenes <- paste("Non-null genes: ", (1 - bhCompare$pNullGenes) * 100, "%", sep = "")
bhCompare$percNonNull <- factor(bhCompare$percNonNull, 
                                levels = paste("Non-zero variants: ", sort(unique(bhCompare$sparsity)) * 100, "%", sep = ""))
bhCompare$Method <- bhCompare$method
bhCompare$Method[bhCompare$method == "poly"] <- "polyhedral"
bhCompare$Method[bhCompare$method == "bh"] <- "BH"
bhCompare$Method[bhCompare$method == "bh"] <- "BH"
ggplot(subset(bhCompare, method != "naive" & test != "skat" & method != "by"), 
       aes(x = R, y = (power), col = Method, linetype = Method, shape = Method)) + 
  geom_point() + geom_line() +
  geom_segment(aes(y = power - 2 * powerSD, yend = power + 2 * powerSD, xend = R)) +
  theme_bw() +
  facet_grid(percGenes ~ percNonNull) + 
  xlab("R_squared") + ylab("Power")
# ggsave(file = "figures/5000genesBHpower.pdf", width = 10, height = 3.5)

# Overall FDR ----------
ggplot(subset(bhCompare, test != "wald" & method != "by"), 
       aes(x = R, y = fdr, col = method, linetype = method, shape = method)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid(sparsity ~ pNullGenes, labeller = "label_both", scales = "free_x") + 
  geom_segment(aes(xend = R, 
                   yend = pmin(fdr + 2 * fdrSD, 1), y = pmax(fdr - 2 * fdrSD, 0))) + 
  geom_hline(yintercept = 0.05) +
  ggtitle("Overall FDR")
# ggsave(file = "figures/5000genesOverallFDR.pdf", width = 10, height = 4)

ggplot(subset(bhCompare, test != "wald" & method != "by"), 
       aes(x = R, y = fdr, col = Method, linetype = Method, shape = Method)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid(percGenes ~ percNonNull) + 
  geom_segment(aes(xend = R, 
                   yend = pmin(fdr + 2 * fdrSD, 1), y = pmax(fdr - 2 * fdrSD, 0))) + 
  geom_hline(yintercept = 0.05)  + 
  xlab("R_squared") + ylab("Overall FDR")
# ggsave(file = "figures/5000genesOverallFDR_A.pdf", width = 10, height = 4)

# Power ----------
power <- do.call("rbind", lapply(results, function(x) x$power))
config <- do.call("rbind", lapply(results, function(x) x$dat))
nGenes <- group_by(config, R, signal, sparsity, pNullGenes, test) %>%
  summarize(discoveries = length(seed) / length(unique(seed)))
config$test <- NULL
power <- cbind(config, power)
power <- gather(power, key = "method", value = "power", 
                hybrid, polyhedral, bbAvgPower)
power <- subset(power, !is.nan(power))
power <- group_by(power, R, sparsity, pNullGenes, method, test, id) %>% 
  summarize(power = mean(power)) %>% 
  group_by(R, sparsity, pNullGenes, method, test) %>% 
  summarize(powerSD = sd(power, na.rm = TRUE) / sqrt(sum(!is.na(power))), 
            power = mean(power, na.rm = TRUE))

ggplot(subset(power, method != "by" & method != "naive" & test != "skat"),
       aes(x = R, y = power, col = method, linetype = method, shape = method)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid(sparsity ~ pNullGenes, labeller = "label_both", scales = "free_x") + 
  geom_segment(aes(xend = R, 
                   yend = pmin(power + 2 * powerSD, 1), y = pmax(power - 2 * powerSD, 0))) + 
  ggtitle("Power of AVG FDR (Over Selected)")
# ggsave(file = "figures/5000genesPower_overSelected.pdf", width = 10, height = 5)

# FDR ----------
fdr <- do.call("rbind", lapply(results, function(x) x$fdr))
fdr <- gather(fdr, key = "method", value = "fdr", hybrid, poly, bb, naive)
fdr <- subset(fdr, !is.nan(fdr))
fdr <- group_by(fdr, R, sparsity, pNullGenes, method, seed, level, test) %>% 
  summarize(fdr = mean(fdr)) %>% 
  group_by(R, sparsity, pNullGenes, method, level, test) %>% 
  summarize(fdrSD = sd(fdr, na.rm = TRUE) / sqrt(sum(!is.na(fdr))), 
            fdr = mean(fdr, na.rm = TRUE))
fdr <- fdr[order(fdr$R), ]
# fdr <- subset(fdr, !(log2(totSNR) == -2 & pNullGenes == 0.975))

normQuant <- qnorm(1 - 0.05 / nrow(fdr))
ggplot(subset(fdr, method != "by" & test != "skat"),
       aes(x = level, y = fdr, col = method, linetype = method, shape = method)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid(R ~ pNullGenes + sparsity, labeller = "label_both") + 
  geom_segment(aes(xend = level, 
                   yend = pmin(fdr + normQuant * fdrSD, 1), 
                   y = pmax(fdr - normQuant * fdrSD, 0))) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  ggtitle("Average FDR")
# ggsave(file = "figures/5000genesfdr.pdf", width = 7, height = 6)

fdr <- subset(fdr, level == 0.05)
fdr <- subset(fdr, test != "skat")
normQuant <- qnorm(1 - 0.05 / nrow(fdr))
fdr$Method <- fdr$method
fdr$Method[fdr$method == "poly"] <- "polyhedral"
fdr$Method[fdr$method == "bb"] <- "BB"
fdr$percNonNull <- paste("Non-zero variants: ", fdr$sparsity * 100, "%", sep = "")
fdr$percGenes <- paste("Non-null genes: ", (1 - fdr$pNullGenes) * 100, "%", sep = "")
fdr$percNonNull <- factor(fdr$percNonNull, 
                                levels = paste("Non-zero variants: ", sort(unique(fdr$sparsity)) * 100, "%", sep = ""))
ggplot(fdr, aes(x = R, y = fdr, col = Method, linetype = Method, shape = Method)) + 
  geom_line() + geom_point() + 
  theme_bw() + 
  facet_grid(percGenes ~ percNonNull) +
  geom_hline(yintercept = 0.05) + 
  geom_segment(aes(xend = R, y = fdr + fdrSD * normQuant, yend = fdr - fdrSD * normQuant)) + 
  xlab("R_squared") + ylab("Average FDR")
# ggsave(file = "figures/5000genes_AverageFDR_A.pdf", width = 10, height = 3.5)

# Power ----------
config <- do.call("rbind", lapply(results, function(x) x$dat))
nGenes <- group_by(config, R, sparsity, pNullGenes, test, id) %>% 
  summarize(discoveries = length(id) / length(unique(id))) %>% 
  group_by(R, sparsity, pNullGenes, test) %>% 
  summarize(discoveries = length(test) / length(unique(test)))
ggplot(nGenes, aes(x = R, y = discoveries, col = test)) + 
  geom_point() + geom_line() +
  facet_grid(sparsity ~ pNullGenes, labeller = "label_both") + 
  theme_bw()
  
power <- do.call("rbind", lapply(results, function(x) x$totpower))
config <- do.call("rbind", lapply(results, function(x) x$dat[1, ]))
config$test <- NULL
power <- rbind(cbind(config, subset(power, test == "wald")), 
                     cbind(config, subset(power, test == "skat")))
power <- gather(power, key = "method", value = "power", hybrid, polyhedral, bbAvgPower)
power <- subset(power, !is.nan(power))
# power <- subset(power, !(log2(totSNR) == -2 & pNullGenes == 0.975))
power <- group_by(power, R, sparsity, pNullGenes, method, test) %>% 
  summarize(powerSD = sd(power) / sqrt(length(power)), power = mean(power))
  # group_by(R, sparsity, pNullGenes, method, test) %>% 
  # summarize(powerSD = sd(power, na.rm = TRUE) / sqrt(sum(!is.na(power))), 
  #           power = mean(power, na.rm = TRUE))

power <- subset(power, method != "by")
power <- subset(power, test != "skat")
power$Method <- power$method
power$Method[power$method == "bbAvgPower"] <- "BB"
power$percNonNull <- paste("Non-zero variants: ", power$sparsity * 100, "%", sep = "")
power$percGenes <- paste("Non-null genes: ", (1 - power$pNullGenes) * 100, "%", sep = "")
power$percNonNull <- factor(power$percNonNull, 
                                levels = paste("Non-zero variants: ", sort(unique(power$sparsity)) * 100, "%", sep = ""))
ggplot(power, aes(x = R, y = power, col = Method, linetype = Method, shape = Method)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid(percGenes ~ percNonNull) + 
  geom_segment(aes(xend = R, 
                   yend = pmin(power + 2 * powerSD, 1), y = pmax(power - 2 * powerSD, 0))) +
  xlab("R_squared") + ylab("Power")
  # ggtitle("Power of AVG FDR (Overall)")
# ggsave(file = "figures/5000genesPower_overall.pdf", width = 10, height = 3.5)
