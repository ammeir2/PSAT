#' Plot results for glmQuadratic Objects
#'
#' @description An interface for plotting the results of a post-selection analysis
#' following an aggregate test. At the moment, methods for plotting
#' estimates and confidence intervals, p-values and stochastic gradietn solution
#' paths are implemented.
#'
#' @param object and object of class \code{mvnQuadratic} or \code{glmQuadratic}.
#'
#' @param type the type of plot to be generated.
#'
#' @param ... additional arguments to be passed to plotting function.
#'
#' @details This functions serves as an interface for plotting the results
#' of post aggregatet testing inference. The function calls the \code{\link{plotEstimates}}
#' function for \code{type} "estimates", the \code{\link{plotSolutionPath}} method for
#' \code{type} "solutionPath" and \code{\link{plotPvalues}} for \code{type} "p-values".
#' See documentation for these functions for details regarding additional optional parameters.
#'
#' @return a \code{\link[ggplot2]{ggplot2}} figure.
#'
#' @seealso \code{\link{mvnQuadratic}}, \code{\link{glmQuadratic}},
#' \code{\link{plotEstimates}}, \code{\link{plotSolutionPath}},
#'  \code{\link{plotPvalues}}
plot.glmQuadratic <- function(object, type = c("estimates", "solution_path", "p-values"), ... ) {
  object$naiveMu <- object$naiveBeta[-1]
  if(!is.null(object$mleBeta)) {
    object$mleMu <- object$mleBeta[-1]
  }

  class(object) <- "mvnQuadratic"
  return(plot(object, type = type, ...))
}

#' Plot results for mvnQuadratic Objects
#'
#' @description see \code{\link{plot.glmQuadratic}} for details.
plot.mvnQuadratic <- function(object, type = c("estimates", "solution-path", "p-values"), ...) {
  type <- type[1]
  if(type == "estimates") {
    return(plotEstimates(object, ...))
  }

  if(type == "solution-path") {
    return(plotSolutionPath(object))
  }

  if(type == "p-values") {
    return(plotPvalues(object, ...))
  }


  stop("Unknown plotting method!")
}

#' Plotting Functions for Post Aggregate Testing Inference Results
#'
#' @description Functions for plotting the results of post aggregate testing
#' analyses.
#'
#' @param true the true values of the estimated means or regression coefficients.
#' Mostly useful for experiments with simulated data.
#'
#' @param threshold a vector of thresholds to be plotted along with the p-values.
#'
#' @param adjust method for adjusting the p-values. See \code{\link[stats]{p.adjust}}
#' for details.
#'
#' @return A \code{\link[ggplot2]{ggplot2}} figure.
#'
#' @describeIn plotEstimates plots estimates and confidence intervals.
plotEstimates <- function(object, true = NULL) {
  naive <- object$naiveMu
  p <- length(naive)
  mle <- object$mleMu
  cis <- object$ci_type
  offset <- 0
  naive <- data.frame(variable = 1:p, estimate = naive,
                      lci = object$naiveCI[, 1], uci = object$naiveCI[, 2],
                      estimate_type = "naive", ci_type = "naive",
                      offset = offset)
  offset <- offset + 1
  if(!is.null(mle)) {
    mle <- data.frame(variable = 1:p, estimate = mle,
                      lci = object$ci[, 1], uci = object$ci[, 2],
                      estimate_type = "mle", ci_type = object$ci_type[1],
                      offset = offset, stringsAsFactors = FALSE)
    offset <- offset + 1
  }

  plotlist <- list()
  slot <- 1
  if(length(object$ci_type) > 1) {
    for(i in 2:length(object$ci_type)) {
      type <- object$ci_type[i]
      if(type == "naive") next
      ci <- getCI(object, type = type)
      plotlist[[slot]] <- data.frame(variable = 1:p, estimate = rep(NA, p),
                                     lci = ci[, 1], uci = ci[, 2],
                                     estimate_type = rep(NA, p), ci_type = type,
                                     offset = rep(offset, p), stringsAsFactors = FALSE)
      slot <- slot + 1
      offset <- offset + 1
    }
  }

  plotdat <- rbind(naive, mle, do.call("rbind", plotlist))

  if(!is.null(true)) {
    if(length(true) != length(object$naiveMu)) {
      stop("Incorrect dimension for underlying true mean!")
    }

    truedat <- data.frame(variable = 1:p, estimate = true,
                          lci = NA, uci = NA, estimate_type = "truth",
                          ci_type = NA, offset = offset,
                          stringsAsFactors = FALSE)
    plotdat <- rbind(plotdat, truedat)
  }

  plotdat$variable <- plotdat$variable + plotdat$offset * 0.12
  ggplot(subset(plotdat, !is.na(ci_type))) +
    geom_segment(aes(x = variable, xend = variable,
                     y = lci, yend = uci, col = ci_type, linetype = ci_type)) +
    theme_bw() + geom_hline(yintercept = 0) +
    geom_point(data = subset(plotdat, !is.na(estimate_type)),
               aes(x = variable, y = estimate, shape = estimate_type)) +
    scale_color_discrete(name = "CI Type") +
    scale_shape_discrete(name = "Estimate Type") +
    scale_linetype_discrete(name = "CI Type") +
    ylab("Estimates / CIs") + xlab("Variable")
}

#' @describeIn plotEstimates plots the solution path of a stochastic
#' gradient method.
plotSolutionPath <- function(object) {
  path <- object$solutionPath
  if(is.null(path)) {
    stop("No solution path, was stochastic optimization performed?")
  }

  path <- reshape2::melt(path)
  names(path) <- c("Iteration", "Variable", "Estimate")
  ggplot(path) + geom_line(aes(x = Iteration, y = Estimate,
                               col = factor(Variable), linetype = factor(Variable))) +
    theme_bw() + scale_linetype_discrete(guide = FALSE) +
    scale_color_discrete(guide = FALSE) +
    geom_hline(yintercept = 0)
}

#' @describeIn plotEstimates plots the selection adjusted p-values.
plotPvalues <- function(object, threshold = 0.1, adjust = "bonferroni") {
  if(is.null(adjust)) {
    adjust <- "none"
  }
  adjust <- adjust[1]
  if(!(adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))) {
    stop("Unknown p-value adjustment method. See documentation for `p.adjust` for details.")
  }

  varnames <- colnames(object$sigma)
  if(is.null(varnames)) {
    varnames <- 1:ncol(object$sigma)
  }
  naive <- data.frame(variable = varnames, p_value = object$naivePval,
                      pvalue_type = "naive", stringsAsFactors = FALSE)

  plotlist <- list()
  slot <- 1
  for(i in 1:length(object$pvalue_type)) {
    type <- object$pvalue_type[i]
    if(type == "naive") next
    pval <- getPval(object, type = type)
    plotlist[[slot]] <- data.frame(variable = varnames, p_value = pval,
                                   pvalue_type = type, stringsAsFactors = FALSE)
    slot <- slot + 1
  }

  plotdat <- rbind(naive, do.call("rbind", plotlist))
  plotdat$variable <- factor(plotdat$variable, levels = varnames)
  plotdat$p_value <- p.adjust(plotdat$p_value, method = adjust)
  plotdat$p_value <- -log10(plotdat$p_value)
  if(adjust == "none") {
    ylab <- c("-log10(p-values)")
  } else {
    ylab <- paste("-log10(", adjust, " adjusted p-values)", sep = "")
  }

  figure <- ggplot(plotdat) +
    geom_point(aes(x = variable, y = p_value, col = pvalue_type,
                                   shape = pvalue_type)) +
    theme_bw() + ylab(ylab)

  if(length(threshold > 0)){
    interceptDat <- data.frame(intercept = -log10(threshold),
                               threshold = factor(threshold))
    figure <- figure +
      geom_hline(data = interceptDat, aes(yintercept = intercept, linetype = threshold))
  }

  return(figure)
}
