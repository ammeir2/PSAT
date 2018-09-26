library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyr)

filenames <- as.list(dir(path = 'simulations/results', pattern="covByRef_H_*"))
filenames <- lapply(filenames, function(x) paste0('simulations/results/', x))
filenames <- as.character(filenames)
results <- vector(length(filenames), mode = "list")
level <- 0.05
pb <- txtProgressBar(min = 0, max = length(filenames), style = 3)
for(i in 1:length(filenames)) {
  res <- readRDS(file = filenames[i])
  sublist <- vector(length(res), mode = "list")
  for(j in 1:length(res)) {
    setting <- res[[j]]
    config <- setting[[1]]
    settingResults <- data
    subsublist <- vector(length(setting) - 1, mode = "list")
    for(k in 2:length(setting)) {
      x <- data.frame(setting[[k]])
      hybridQvals <- p.adjust(x$hybrid, method = "BH")
      polyQvals <- p.adjust(x$polyhedral, method = "BH")
      naiveQvals <- p.adjust(x$naive, method = "BH")
      nonzero <- x$true != 0
      hybridFDR <- sum(hybridQvals < level & !nonzero) / max(sum(hybridQvals < level), 1)
      polyFDR <- sum(polyQvals < level & !nonzero) / max(sum(polyQvals < level), 1)
      naiveFDR <- sum(naiveQvals < level & !nonzero) / max(sum(naiveQvals < level), 1)
      fdr <- data.frame(hybrid = hybridFDR, polyhedral = polyFDR, naive = naiveFDR, type = "FDR")
      hybridPower <- sum(hybridQvals < level & nonzero) / sum(nonzero)
      polyPower <- sum(polyQvals < level & nonzero) / sum(nonzero)
      naivePower <- sum(naiveQvals < level & nonzero) / sum(nonzero)
      power <- data.frame(hybrid = hybridPower, polyhedral = polyPower, naive = naivePower, type = "power")
      subsublist[[k - 1]] <- rbind(fdr, power)
    }
    sublist[[j]] <- do.call("rbind", subsublist) %>% cbind(t(config), .)
  }
  results[[i]] <- do.call("rbind", sublist)
  setTxtProgressBar(pb, i)
}
close(pb)
dims <- sapply(results, dim)
# results <- results[dims[2, ] == 13]
results <- do.call("rbind", results)
# saveRDS(results, file = "simulations/results/covByRef_processed.rds")

# False Detection Rate ----- 
fdr <- subset(results, type == "FDR") %>% dplyr::select(-type) %>% 
  gather(key = method, value = fdr, hybrid, polyhedral, naive) %>% 
  group_by(n, refSize, p, snr, sparsity, rho, pvalThreshold, method) %>%
  summarize(R_squared = mean(R_squared), fdrSD = sd(fdr) / sqrt(length(fdr)), fdr = mean(fdr)) 
nq <- qnorm(1 - 0.05 / nrow(fdr))
fdr$R_squared <- round(fdr$R_squared, 4)
ggplot(fdr, aes(x = factor(R_squared), y = fdr, col = method, linetype = method, shape = method)) +
  geom_point() + 
  geom_hline(yintercept = level) + 
  theme_bw() +
  geom_segment(aes(xend = factor(R_squared), y = pmax(fdr - nq * fdrSD, 0), yend = pmin(fdr + nq * fdrSD, 1))) +
  facet_grid(p + refSize ~ rho + sparsity, labeller = "label_both", scales = "free") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file = "figures/covByRef_FDR_A.pdf", width = 12, height = 8)

# Power ---------
power <- subset(results, type == "power") %>% dplyr::select(-type) %>% 
  gather(key = method, value = power, hybrid, polyhedral, naive) %>% 
  group_by(refSize, p, snr, sparsity, rho, pvalThreshold, method) %>% 
  summarize(powerSD = sd(power) / sqrt(length(power)), power = mean(power), 
            R_squared = mean(R_squared)) 
nq <- qnorm(1 - 0.05 / nrow(power))
power$R_squared <- round(power$R_squared, 4)
ggplot(power, aes(x = factor(R_squared), y = power, col = method, linetype = method, shape = method)) +
  geom_point() + 
  theme_bw() +
  geom_segment(aes(xend = factor(R_squared), y = pmax(power - nq * powerSD, 0), yend = pmin(power + nq * powerSD, 1))) +
  facet_grid(p + refSize ~ rho + sparsity, labeller = "label_both", scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file = "figures/covByRef_power_A.pdf", width = 12, height = 8)

