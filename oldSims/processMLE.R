library(magrittr)
library(ggplot2)
library(dplyr)
library(reshape2)

filenames <- as.list(dir(path = 'oldSims/results', pattern="mle_sim_A_*"))
filenames <- lapply(filenames, function(x) paste0('oldSims/results/', x))
filenames <- as.character(filenames)
results <- list(length(filenames))
for(i in 1:length(filenames)) {
  results[[i]] <- readRDS(filenames[[i]])
}

results <- do.call("c", results)
results <- do.call("c", results)
colnames <- names(results[[1]]$params)
results <- sapply(results, function(x) matrix(x[[1]], nrow = 1)) %>% t()
results <- data.frame(results)
names(results) <- colnames
results$m <- NULL
results$selectProb <- NULL
results$chisqPval <- NULL
results$coefType <- NULL
results$noiseType <- NULL
results$MAFthreshold <- NULL
results$pthreshold <- NULL
results$t2type <- NULL

library(reshape2)
library(ggplot2)
library(dplyr)
results <- melt(results, id = c("n", "p", "sparsity", "signal", "R_squared"))
names(results)[6:7] <- c("method", "mse")
results <- summarize(group_by(results, n, p, signal, method),
                     R_squared = round(mean(R_squared), 6),
                     msesd = sd(mse) / length(mse),
                     mse = mean(mse))
results$lci <- results$mse - 2 * results$msesd
results$uci <- results$mse + 2 * results$msesd
results$methodTemp <- results$method
results$method <- as.character(results$method)
results$method[results$method == "mle"] <- "conditional mle"
results$method[results$method == "naive"] <- "naive mle"
results$method <- factor(results$method,
                         levels = c("naive mle", "conditional mle"))
results$Dimension <- results$p
results$snr <- results$signal

# pdf("figures/mleplot_B.pdf",pagecentre=T, width=8,height=2.5 ,paper = "special")
ggplot(subset(results, signal != 0.05),
       aes(x = n / 1000, y = mse,
                    col = method, shape = method, linetype = method)) +
  theme_bw() + facet_grid(R_squared ~ Dimension, labeller = "label_both") +
  geom_point() + geom_line() +
  geom_segment(aes(xend = n / 1000, y = lci, yend = uci)) +
  xlab("Sample Size / 1000") + ylab("RMSE")
# dev.off()





