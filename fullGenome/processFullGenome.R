library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)

filenames <- as.list(dir(path = 'fullGenome/results', pattern="fullGenome_C_*"))
filenames <- lapply(filenames, function(x) paste0('fullGenome/results/', x))
filenames <- as.character(filenames)
results <- vector(length(filenames), mode = "list")
power <- list()
cover <- list()
for(i in 1:length(filenames)) {
  result <- readRDS(filenames[[i]])
  if(length(result[[1]]) == 1) {
    next
  }
  dat <- lapply(result, function(x) x$config) %>% do.call("rbind", .)
  power <- lapply(result, function(x) x$power) %>% do.call("rbind", .)
  cover <- lapply(result, function(x) x$cover) %>% do.call("rbind", .)
  # print(c(i, nrow(dat), nrow(power)))
  if(nrow(dat) != nrow(power)) {
    print("hello")
  }
  results[[i]] <- list(dat = dat, power = power, cover = cover)
}
results <- results[!sapply(results, is.null)]

# Power ----------
power <- do.call("rbind", lapply(results, function(x) x$power))
config <- do.call("rbind", lapply(results, function(x) x$dat))
power <- cbind(config, power)
nGenes <- group_by(config, totSNR, signal, sparsity, pNullGenes) %>% 
  summarize(discoveries = length(seed) / length(unique(seed)))
power <- melt(power, id = colnames(power)[-(13:15)])
names(power)[13:14] <- c("method", "power")
power <- subset(power, !is.nan(power))
power <- group_by(power, totSNR, sparsity, pNullGenes, method, seed) %>% 
  summarize(power = mean(power)) %>% 
  group_by(totSNR, sparsity, pNullGenes, method) %>% 
  summarize(powerSD = sd(power, na.rm = TRUE) / sqrt(sum(!is.na(power))), 
            power = mean(power, na.rm = TRUE))

temp <- subset(power, totSNR != 0.125)
temp <- subset(temp, !(totSNR == 0.25 & pNullGenes == 0.975))
ggplot(temp, aes(x = log2(totSNR), y = power, col = method, linetype = method, shape = method)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid(sparsity ~ pNullGenes, labeller = "label_both", scales = "free_x") + 
  geom_segment(aes(xend = log2(totSNR), 
                   yend = pmin(power + 2 * powerSD, 1), y = pmax(power - 2 * powerSD, 0)))
# ggsave(file = "figures/5000genesPower.pdf", width = 5, height = 2.5)

# Cover ----------
cover <- do.call("rbind", lapply(results, function(x) x$cover))
config <- do.call("rbind", lapply(results, function(x) x$dat))
cover <- cbind(config, cover)
cover <- melt(cover, id = colnames(cover)[-(13:15)])
names(cover)[13:14] <- c("method", "cover")
cover <- subset(cover, !is.nan(cover))
cover <- group_by(cover, totSNR, sparsity, pNullGenes, method, seed) %>% 
  summarize(cover = mean(cover)) %>% 
  group_by(totSNR, sparsity, pNullGenes, method) %>% 
  summarize(coverSD = sd(cover, na.rm = TRUE) / sqrt(sum(!is.na(cover))), 
            cover = mean(cover, na.rm = TRUE))

temp <- subset(cover, totSNR != 0.125)
temp <- subset(temp, !(totSNR == 0.25 & pNullGenes == 0.975))
ggplot(temp, aes(x = totSNR, y = cover, col = method, linetype = method, shape = method)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_grid(sparsity ~ pNullGenes, labeller = "label_both") + 
  geom_segment(aes(xend = totSNR, 
                   yend = pmin(cover + 2 * coverSD, 1), 
                   y = pmax(cover - 2 * coverSD, 0))) + 
  geom_hline(yintercept = 0.95)
# ggsave(file = "figures/5000genesCover.pdf", width = 5, height = 2.5)






