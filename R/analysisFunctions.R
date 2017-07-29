getCI <- function(object, type = NULL, confidence_level = NULL, switchTune = NULL) {
  # Checking arguments ---------
  if(!(class(object) %in% c("mvnQuadratic", "glmQuadratic"))) {
    stop("Invalid object type!")
  }

  if(class(object) == "glmQuadratic") {
    object$muhat <- object$betahat[-1]
    object$naiveMu <- object$naiveBeta[-1]
  }

  if(is.null(type)) {
    type <- object$ci_type[1]
  } else {
    type <- type[1]
  }

  if(is.null(switchTune)) {
    switchTune <- object$switchTune
  }

  if(is.null(confidence_level)) {
    confidence_level <- object$confidence_level
    recompute <- FALSE
  } else if(confidence_level == object$confidence_level) {
    recompute <- FALSE
  } else {
    recompute <- TRUE
  }

  # Polyhedral ----------
  if(type == "polyhedral") {
    if(recompute | is.null(object$polyCI)) {
      return(getPolyCI(object$naiveMu, object$sigma, object$testMat,
                       object$threshold, confidence_level)$ci)
    } else {
      return(object$polyCI)
    }
  }

  # Regime switching --------------
  if(type == "switch") {
    if(switchTune != object$switchTune) {
      recompute <- TRUE
    }
    if(is.null(object$nullDist)) {
      stop("Please re-run fitting function with `switch` or `global-null` ci_type options.")
    } else if(!recompute) {
      return(object$switchCI)
    } else {
      if(confidence_level != object$confidence_level) {
        naiveCI <- getNaiveCI(object$naiveMu, object$sigma, confidence_level)
      } else {
        naiveCI <- object$naiveCI
      }
      switch <- getSwitchCI(object$naiveMu, object$muhat, object$nullDist,
                            object$testMat, switchTune, confidence_level,
                            object$pthreshold, object$quadlam, naiveCI)
      return(switch)
    }
  }

  # Naive -------------
  if(type == "naive") {
    if(!recompute) {
      return(object$naiveCI)
    } else {
      return(getNaiveCI(object$naiveMu, object$sigma, confidence_level))
    }
  }

  # Global Null -------------
  if(type == "global-null") {
    if(is.null(object$nullDist)) {
      stop("Please re-run fitting function with `switch` or `global-null` ci_type options.")
    } else if(!recompute) {
      return(object$nullCI)
    } else {
      return(getNullCI(object$muhat, object$nullDist, confidence_level))
    }
  }

  # MLE CI ----------------------
  if(type == "mle") {
    if(is.null(object$mleDist)) {
      stop("Please re-run fitting function with `mle` ci_type option.")
    } else if(!recompute) {
      return(object$mleCI)
    } else {
      return(getMleCI(object$muhat, object$mleDist, confidence_level))
    }
  }

  # Oopps --------
  stop("Unsupported CI type!")
}

getPval <- function(object, type = NULL) {
  # Checking arguments ---------
  if(!(class(object) %in% c("mvnQuadratic", "glmQuadratic"))) {
    stop("Invalid object type!")
  }

  if(class(object) == "glmQuadratic") {
    object$muhat <- object$betahat[-1]
    object$naiveMu <- object$naiveBeta[-1]
  }

  type <- type[1]
  if(is.null(type)) {
    type <- object$pvalue_type[1]
  }

  # Polyhedral ----------------------
  if(type == "polyhedral") {
    if(is.null(object$polyPval)) {
      return(getPolyCI(object$naiveMu, object$sigma, object$testMat,
                object$threshold, confidence_level, FALSE)$pval)
    } else {
      return(object$polyPval)
    }
  }

  # Hybrid ------------------
  if(type == "hybrid") {
    if(is.null(object$nullSample)) {
      stop("Please re-run fitting function with `hybrid` or `global-null` pvalue_type option.")
    }
    if(is.null(object$polyPval)) {
      polyPval <- getPolyCI(object$naiveMu, object$sigma, object$testMat,
                            object$threshold, confidence_level, FALSE)$pval
    } else {
      polyPval <- object$polyPval
    }
    return(pmin(2 * pmin(polyPval, object$nullPval), 1))
  }

  # Global Null ----------------
  if(type == "global-null") {
    if(!is.null(object$nullPval)) {
      return(object$nullPval)
    } else {
      stop("Please re-run fitting function with `hybrid` or `global-null` pvalue_type option.")
    }
  }

  # Naive ---------------
  if(type == "naive") {
    return(object$naivePval)
  }

  # MLE --------------
  if(type == "mle") {
    if(is.null(object$mlePval)) {
      stop("Please re-run fitting function with `mle` pvalue_type option.")
    } else {
      return(object$mlePval)
    }
  }

  # Oopps --------
  stop("Unsupported pvalue type!")
}

coef.mvnQuadratic <- function(object, type = NULL) {
  if(is.null(type)) {
    return(object$muhat)
  }

  if(type == "mle") {
    if(is.null(object$mleMu)) {
      stop("Please re-run fitting routine with `mle' estimate_type option.")
    } else {
      return(object$mleMu)
    }
  }

  if(type == "naive") {
    return(object$naiveMu)
  }

  stop("Unsupported estimate type!")
}

coef.glmQuadratic <- function(object, type = NULL) {
  if(is.null(type)) {
    return(object$betahat)
  }

  if(type == "mle") {
    if(is.null(object$mleBeta)) {
      stop("Please re-run fitting routine with `mle' estimate_type option.")
    } else {
      return(object$mleBeta)
    }
  }

  if(type == "naive") {
    return(object$naiveBeta)
  }

  stop("Unsupported estimate type!")
}

predict.mvnQuadratic <- function(object) {
  return(object$muhat)
}

predict.glmQuadratic <- function(object, newX = NULL) {
  if(is.null(newX)) {
    X <- object$X
  } else {
    X <- newX
  }

  X <- cbind(1, X)
  return(as.numeric(X %*% object$betahat))
}

plot.glmQuadratic <- function(object, type = c("estimates", "solution_path"), ... ) {
  object$naiveMu <- object$naiveBeta[-1]
  if(!is.null(object$mleBeta)) {
    object$mleMu <- object$mleBeta[-1]
  }

  class(object) <- "mvnQuadratic"
  return(plot(object, type = type, ...))
}

plot.mvnQuadratic <- function(object, type = c("estimates", "solution-path"), ...) {
  require(ggplot2)
  type <- type[1]
  if(type == "estimates") {
    return(plotEstimates(object, ...))
  }

  if(type == "solution-path") {
    return(plotSolutionPath(object))
  }

  stop("Unknown plotting method!")
}

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
  }

  plotdat <- rbind(plotdat, truedat)
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

summary.glmQuadratic <- function(object, estimate_type = NULL, pvalue_type = NULL,
                                 ci_type = NULL, confidence_level = NULL) {
  est <- coef(object, type = estimate_type)
  if(is.null(estimate_type)) {
    estimate_type <- object$estimate_type[1]
  }
  if(estimate_type == "mle") {
    estimate_type <- "MLE"
  } else {
    estimate_type <- "Naive"
  }

  pvalue <- getPval(object, type = pvalue_type)
  if(is.null(pvalue_type)) {
    pvalue_type <- object$pvalue_type[1]
  }
  if(pvalue_type == "naive") {
    pvalue_type <- "Naive"
  } else if(pvalue_type == "polyhedral") {
    pvalue_type <- "Polyhedral"
  } else if(pvalue_type == "hybrid") {
    pvalue_type <- "Hybrid"
  } else if(pvalue_type == "global-null") {
    pvalue_type <- "Global Null"
  } else {
    pvalue_type <- "MLE"
  }

  if(is.null(confidence_level)) {
    confidence_level <- object$confidence_level
  }
  ci <- getCI(object, type = ci_type, confidence_level = confidence_level)
  if(is.null(ci_type)) {
    ci_type <- object$ci_type[1]
  }
  if(ci_type == "naive") {
    ci_type <- "Naive"
  } else if(ci_type == "switch") {
    ci_type <- "Regime Switching"
  } else if(ci_type == "global-null"){
    ci_type <- "Global Null"
  } else if(ci_type == "polyhedral") {
    ci_type <- "Polyhedral"
  } else {
    ci_type <- "MLE"
  }

  intsd <- object$interceptsd
  intpval <- 2 * pnorm(-abs(est[1] / intsd))
  intci <- est[1] + c(-1, 1) * intsd * qnorm(1 - confidence_level / 2)

  coefficients <- data.frame(estimates = est, lCI = c(intci[1], ci[, 1]),
                             uCI = c(intci[2], ci[, 2]),
                             pvalue = c(intpval, pvalue))
  rownames <- colnames(object$X)
  if(is.null(rownames)) {
    rownames <- paste("V", 0:(length(est) - 1), sep ="")
    rownames[1] <- "intercept"
  }
  rownames(coefficients) <- rownames

  sum <- list()
  sum$coefficients <- coefficients
  sum$confidence_level <- confidence_level
  sum$ci_type <- ci_type
  sum$pvalue_type <- pvalue_type
  sum$estimate_type <- estimate_type
  class(sum) <- "glmQuadratic_summary"
  return(sum)
}

print.glmQuadratic_summary <- function(sum) {
  cat("Results for Post-Aggregate Testing Analysis \n \n")
  cat("Estimation Method: ", sum$estimate_type, "\n")
  cat("P-value Type: ", sum$pvalue_type, "\n")
  cat("Confidence Interval Type: ", sum$ci_type, "      Confidence Level:", sum$confidence_level, "\n\n")
  cat("Table of Coefficients:\n")
  coefs <- sum$coefficients
  coefs[, 1:3] <- round(coefs[, 1:3], 5)
  names(coefs) <- c("Estimate", "Lower-CI", "Upper-CI", "P-Value")

  p <- nrow(coefs)
  siglevel <- rep("", p)
  siglevel[coefs[, 4] < 0.1] <- "."
  siglevel[coefs[, 4] < 0.05] <- "*"
  siglevel[coefs[, 4] < 0.01] <- "**"
  siglevel[coefs[, 4] < 0.001] <- "***"
  coefs <- cbind(coefs, siglevel)
  names(coefs)[5] <- ""

  print(coefs)
  cat("\n Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
}






