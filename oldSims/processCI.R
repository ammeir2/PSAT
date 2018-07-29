library(magrittr)
library(ggplot2)
library(dplyr)
library(reshape2)

filenames <- as.list(dir(path = 'oldSims/results', pattern="ci_sim_A_*"))
filenames <- lapply(filenames, function(x) paste0('oldSims/results/', x))
filenames <- as.character(filenames)
results <- list(length(filenames))
for(i in 1:length(filenames)) {
  results[[i]] <- readRDS(filenames[[i]])
}

results <- do.call("c", results)
simResults <- results

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
                      R_squared = iterResult[[1]]$params[["R_squared"]],
                      cover = cover, measure = "cover")
  det <- data.frame(p = iterResult[[1]]$params[3],
                    signal = iterResult[[1]]$params[[5]],
                    t2 = iterResult[[1]]$params[[12]],
                    sparisty = iterResult[[1]]$params[[4]],
                    R_squared = iterResult[[1]]$params[["R_squared"]],
                    det = det, measure = "signDet")

  names(cover) <- c("p", "signal", "t2", "sparsity", "R_squared", "naive", "switch", "poly", "measure")
  names(det) <- c("p", "signal", "t2", "sparsity", "R_squared", "naive", "switch", "poly", "measure")
  coverList[[i]] <- cover
  detList[[i]] <- det
}

library(dplyr)
library(reshape2)
library(ggplot2)
results <- rbind(do.call("rbind", coverList), do.call("rbind", detList))
results <- melt(results, id = c("p", "signal", "t2", "sparsity", "measure", "R_squared"))
names(results)[7:8] <- c("method", "error")
switch <- subset(results, method == "switch")
switch$method <- NULL
switch$method[switch$t2 == 1] <- "switch-sqrd"
switch$method[switch$t2 == 2] <- "switch-half"
switch$t2 <- NULL
results <- subset(results, method != "switch")
results$t2 <- NULL
results <- rbind(results, switch)
results <- summarize(group_by(results, p, signal, sparsity, method, measure),
                     R_squared = mean(R_squared),
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

# pdf("figures/ciplotWpower_B.pdf",pagecentre=T, width=8, height=3.5 ,paper = "special")
results$R_squared <- log2(results$R_squared)
ggplot(subset(results, sparsity != 10), aes(x = R_squared, y = error, col = method, linetype = method)) +
  geom_line() + geom_point(aes(shape = method)) +
  facet_grid(measureFacet ~ sparseFacet, scales = "free") +
  geom_segment(aes(x = R_squared, xend = R_squared, y = lci, yend = uci)) +
  theme_bw() + geom_hline(aes(yintercept = intercept)) +
  #geom_hline(aes(yintercept = 1 - (1 - intercept) / 2), linetype = 2) +
  ylab("") + xlab("log2(R_squared)")
# dev.off()
