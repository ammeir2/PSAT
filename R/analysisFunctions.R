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






