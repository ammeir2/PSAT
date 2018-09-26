library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyr)

filenames <- as.list(dir(path = 'simulations/results', pattern="unknownVar_FIXED_A_*"))
filenames <- as.list(dir(path = 'simulations/results', pattern="unknownVar_FIXED_B_*"))
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
      pvals <- data.frame(hybrid = x$hybrid, polyhedral = x$polyhedral, naive = x$naive, 
                          knownHybrid = x$knownHybrid, knownPolyhedral = x$knownPolyhedral)
                          # nVarHybrid = x$naiveHybrid, nVarPolyhedral = x$naivePolyhedral)
      qvals <- apply(pvals, 2, p.adjust, method = "BH")
      nonzero <- x$true != 0
      fdr <- apply(qvals, 2, function(x) sum(x < 0.05 & !nonzero) / max(sum(x < 0.05), 1)) %>% 
        matrix(nrow = 1) %>% data.frame() %>%  cbind(type = "FDR")
      names(fdr)[1:ncol(qvals)] <- colnames(qvals)
      power <- apply(qvals, 2, function(x) sum(x < 0.05 & nonzero) / max(sum(nonzero), 1)) %>% 
        matrix(nrow = 1) %>% data.frame() %>% cbind(type = "power")
      names(power)[1:ncol(qvals)] <- colnames(qvals)
      subsublist[[k - 1]] <- rbind(fdr, power)
    }
    sublist[[j]] <- do.call("rbind", subsublist) %>% cbind(t(config), .)
  }
  results[[i]] <- do.call("rbind", sublist)
  setTxtProgressBar(pb, i)
}
close(pb)
results <- do.call("rbind", results)
# saveRDS(results, file = "simulations/results/unknownVar_processed_A.rds")
# results <- readRDS(file = "simulations/results/unknownVar_processed_A.rds")

# False Detection Rate ----- 
fdr <- subset(results, type == "FDR") %>% dplyr::select(-type, -naive) %>% 
  gather(key = method, value = fdr, hybrid, polyhedral, knownHybrid, knownPolyhedral) %>% 
  group_by(n, p, snr, sparsity, rho, pvalThreshold, method) %>%
  summarize(fdrSD = sd(fdr) / sqrt(length(fdr)), fdr = mean(fdr), R_squared = round(mean(R_squared), 2)) 
nq <- qnorm(1 - 0.05 / nrow(fdr))
fdr <- subset(fdr)
ggplot(fdr, aes(x = R_squared, y = fdr, col = method, linetype = method, shape = method)) +
  geom_point() + geom_line() + 
  geom_hline(yintercept = level) + 
  theme_bw() +
  geom_segment(aes(xend = R_squared, y = pmax(fdr - nq * fdrSD, 0), yend = pmin(fdr + nq * fdrSD, 1))) +
  facet_grid(p + n ~ rho + sparsity, labeller = "label_both", scales = "free")
# ggsave(file = "figures/unknownVar_FDR_A.pdf", width = 10, height = 6)

fdr <- subset(fdr, sparsity == 3 & p == 20)
fdr <- subset(fdr, method != "nVarHybrid" & method != "nVarPolyhedral")
fdr <- subset(fdr, rho == 0.5)
nq <- qnorm(1 - 0.05 / nrow(fdr))
fdr$Method <- fdr$method
fdr$Method[fdr$method == "knownHybrid"] <- "hybrid-known"
fdr$Method[fdr$method == "knownPolyhedral"] <- "polyhedral-known"
ggplot(fdr, aes(x = R_squared, y = fdr, col = Method, linetype = Method, shape = Method)) +
  geom_point() + geom_line() + 
  geom_hline(yintercept = level) + 
  theme_bw() +
  geom_segment(aes(xend = R_squared, y = pmax(fdr - nq * fdrSD, 0), yend = pmin(fdr + nq * fdrSD, 1))) +
  facet_grid(rho ~ n, labeller = "label_both") + 
  ylab("FDR") + ylim(0, 0.1)
# ggsave(file = "figures/unknownVar_FDR_A.pdf", width = 10, height = 3)

# Power ---------
power <- subset(results, type == "power") %>% dplyr::select(-type, - naive) %>% 
  gather(key = method, value = power, hybrid, polyhedral, knownHybrid, knownPolyhedral) %>% 
  group_by(n, p, snr, sparsity, rho, pvalThreshold, method) %>% 
  summarize(R_squared = mean(R_squared), powerSD = sd(power) / sqrt(length(power)), power = mean(power)) 
nq <- qnorm(1 - 0.05 / nrow(power))
ggplot(power, aes(x = R_squared, y = power, col = method, linetype = method, shape = method)) +
  geom_point() + geom_line() + 
  theme_bw() +
  geom_segment(aes(xend = R_squared, y = pmax(power - nq * powerSD, 0), yend = pmin(power + nq * powerSD, 1))) +
  facet_grid(p + n ~ rho + sparsity, labeller = "label_both", scales = "free")
# ggsave(file = "figures/unknownVar_power_A.pdf", width = 10, height = 6)

power <- subset(power, sparsity == 3)
power <- subset(power, method != "nVarHybrid" & method != "nVarPolyhedral")
power <- subset(power, rho == 0.5)
nq <- qnorm(1 - 0.05 / nrow(power))
power$Method <- power$method
power$Method[power$method == "knownHybrid"] <- "hybrid-known"
power$Method[power$method == "knownPolyhedral"] <- "polyhedral-known"
ggplot(power, aes(x = R_squared, y = power, col = Method, linetype = Method, shape = Method)) +
  geom_point() + geom_line() + 
  theme_bw() +
  geom_segment(aes(xend = R_squared, y = pmax(power - nq * powerSD, 0), yend = pmin(power + nq * powerSD, 1))) +
  facet_grid(rho ~ n, labeller = "label_both", scales = "free")
# ggsave(file = "figures/unknownVar_power_A.pdf", width = 10, height = 3)

# Joint figure -----------------
fdr$error <- "FDR"
power$error <- "Power"
fdr$value <- fdr$fdr
power$value <- power$power
fdr$sd <- fdr$fdrSD
power$sd <- power$powerSD
fdr$intercept <- 0.05
power$intercept <- NA
forplot <- rbind(fdr, power)
nq <- 2#qnorm(1 - 0.05 / nrow(forplot))
forplot$sampSize <- paste("Sample Size:", forplot$n)
forplot$sampSize <- factor(forplot$sampSize, 
                           levels = paste("Sample Size:", sort(unique(forplot$n))))
ggplot(forplot, aes(x = R_squared, y = value, col = Method, linetype = Method, shape = Method)) + 
  geom_point() + geom_line() +
  geom_segment(aes(xend = R_squared, yend = value - nq * sd, y = value + nq * sd)) + 
  theme_bw() + 
  facet_grid(error ~ sampSize, scales = "free") + 
  ylab("") + 
  geom_hline(aes(yintercept = intercept), linetype = 2)
# ggsave(file = "figures/unknownVar_bothFigure.pdf", width = 10, height = 3.5)

# Numbers
fdr %>% subset(method == "polyhedral") %>% subset(rho == 0) %>% subset(n == 50) %>% subset(snr == 0)
fdr %>% subset(method == "polyhedral") %>% subset(rho == 0) %>% subset(n == 200) %>% subset(snr == 0)
fdr %>% subset(method == "hybrid") %>% subset(rho == 0) %>% subset(n == 50) %>% subset(snr == 0)
fdr %>% subset(method == "hybrid") %>% subset(rho == 0) %>% subset(n == 200) %>% subset(snr == 0)

fdr %>% subset(method == "polyhedral") %>% subset(rho == 0.5) %>% subset(n == 50) %>% subset(snr == 0)
fdr %>% subset(method == "polyhedral") %>% subset(rho == 0.5) %>% subset(n == 200) %>% subset(snr == 0)
fdr %>% subset(method == "hybrid") %>% subset(rho == 0.5) %>% subset(n == 50) %>% subset(snr == 0)
fdr %>% subset(method == "hybrid") %>% subset(rho == 0.5) %>% subset(n == 200) %>% subset(snr == 0)
