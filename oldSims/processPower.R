library(magrittr)
library(ggplot2)
library(dplyr)
library(reshape2)

filenames <- as.list(dir(path = 'oldSims/results', pattern="power_sim_A_*"))
filenames <- lapply(filenames, function(x) paste0('oldSims/results/', x))
filenames <- as.character(filenames)
results <- list(length(filenames))
for(i in 1:length(filenames)) {
  results[[i]] <- readRDS(filenames[[i]])
}

results <- do.call("c", do.call("c", results))
processResults <- function(res, level = 0.05) {
  par <- res$params
  sparsity <- par[[4]]

  poly <- p.adjust(c(res$polyAlt, res$conditionalNullPvals), method = "BH")
  poly <- mean(poly[1:sparsity] < level)
  naive <- p.adjust(c(res$naiveAltPvals, res$naiveNullPvals), method = "BH")
  naive <- mean(naive[1:sparsity] < level)
  null <- p.adjust(c(res$conditionalAltPvals, res$conditionalNullPvals), method = "BH")
  null <- mean(null[1:sparsity] < level)
  hybrid <- p.adjust(c(res$altNotScore, res$nullNotScore), method = "BH")
  hybrid <- mean(hybrid[1:sparsity] < level)
  result <- cbind(t(replicate(4, par)),
                  data.frame(power = c(poly, naive, null, hybrid),
                             method = c("polyhedral", "naive", "global-null", "hybrid")))
  return(result)
}

results <- lapply(results, processResults)
results <- data.table::rbindlist(results)

library(dplyr)
library(ggplot2)
results <- summarize(group_by(results, n, p, sparsity, signal, method),
                     R_squared = mean(R_squared),
                     sdpower = sd(power) / sqrt(length(power)),
                     power = mean(power))
results$lci <- results$power - results$sdpower * 2
results$uci <- results$power + results$sdpower * 2

results$method <- as.character(results$method)
results$method <- factor(results$method, levels = c("naive", "polyhedral", "global-null", "hybrid"))
names(results)[names(results) == "signal"] <- "snr"

# pdf("figures/powerplot_B.pdf",pagecentre=T, width=8,height=2.5 ,paper = "special")
results$R_squared <- round(results$R_squared, 3)
ggplot(subset(results, snr %in% c(0.032, 0.064, 0.128, 0.256)), aes(x = sparsity, y = power, col = method, linetype = method, shape = method)) +
  geom_line() + geom_point() + theme_bw() +
  ylim(0, 1) +
  geom_segment(aes(xend = sparsity, y = lci, yend = uci)) +
  facet_grid(. ~ R_squared, labeller = "label_both") +
  xlab("Sparsity") + ylab("Power")
# dev.off()





