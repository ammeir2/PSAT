#' Retrieve Selection Adjusted Confidence Intervals
#'
#' @description A function for retrieving selection adjusted confidence intervals
#' obtained after aggregate testing using the \code{\link{mvnQuadratic}}
#' or \code{\link{glmQuadratic}} functions.
#'
#' @param object an object of class \code{\link{mvnQuadratic}} or
#' \code{\link{glmQuadratic}}.
#'
#' @param type the type of confidence interval to contruct, see the
#' the documentation of \code{\link{mvnQuadratic}} for details.
#'
#' @param confidence_level the desired confidence level for the confidence
#' intervals
#'
#' @param switchTune tuning parameter for Regime Switching confidence intervals.
#'
#' @details This function retrieves confidence intervals from a post-aggregate testing analysis.
#' if \code{type} is left as \code{NULL} then the confidence interval is the first listed in the
#' \code{ci_type} field of the original function call. If \code{switchTune} is left
#' as \code{NULL} then the method used in the original function call will be used.
#'
#' @return a set of confidence intervals for the normal means/regression
#' coefficients.
getCI <- function(object, type = NULL, confidence_level = NULL, switchTune = NULL) {
  # Checking arguments ---------
  if(!any(class(object) %in% c("mvnQuadratic", "psatGLM", "mvnLinear"))) {
    stop("Invalid object type!")
  }

  aggregateTest <- object$testType

  if("psatGLM" %in% class(object)) {
    object$muhat <- object$betahat[-1]
    object$naiveMu <- object$naiveBeta[-1]
  }

  if(is.null(type)) {
    type <- object$ci_type[1]
  } else {
    type <- type[1]
  }

  if(is.null(confidence_level)) {
    confidence_level <- 1 - object$confidence_level
    recompute <- FALSE
  } else if(confidence_level == object$confidence_level) {
    recompute <- FALSE
    confidence_level <- 1 - confidence_level
  } else {
    recompute <- TRUE
    confidence_level <- 1 - confidence_level
  }

  if(is.null(switchTune)) {
    switchTune <- confidence_level^2
  }

  # Polyhedral ----------
  if(type == "polyhedral") {
    if(recompute | is.null(object$polyCI)) {
      return(getPolyCI(object$naiveMu, object$sigma, object$testMat,
                       object$threshold, confidence_level,
                       test = aggregateTest)$ci)
    } else {
      return(object$polyCI)
    }
  }

  # Regime switching --------------
  if(type == "switch") {
    if(switchTune != object$switchTune) {
      recompute <- TRUE
    }
    if(!recompute & !is.null(object$switchCI)) {
      return(object$switchCI)
    } else {
      if(aggregateTest == "linear") {
        testmat <- object$contrast
      } else {
        testmat <- object$testMat
      }
      switch <- getSwitchCI(object$naiveMu, object$sigma,
                            testmat, object$threshold,
                            object$pthreshold, confidence_level,
                            object$quadlam, t2 = switchTune * object$pthreshold,
                            object$testStat, object$hybridPval,
                            object$trueHybrid, test = aggregateTest)
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
    if(!recompute & !is.null(object$nullCI)) {
      return(object$nullCI)
    } else {
      if(object$nullMethod == "zero-quantile") {
        return(getNullCI(object$muhat, object$nullDist, confidence_level))
      } else {
        if(aggregateTest == "quadratic") {
          ci <- quadraticRB(object$naiveMu, object$sigma,
                            object$testMat, object$threshold,
                            confidence_level, computeFull = TRUE)
        } else if(aggregateTest == "linear") {
          ci <- linearRB(object$naiveMu, object$sigma,
                         object$contrast, object$threshold,
                         object$sigma, computeFull = TRUE)
        }
        return(ci)
      }
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

  if(type == "hybrid") {
    if(!recompute & !is.null(object$hybridCI)) {
      return(object$hybridCI)
    } else {
      ci <- getHybridCI(object$naiveMu, object$sigma,
                        object$testMat, object$threshold,
                        object$pthreshold, confidence_level,
                        object$hybridPval, object$trueHybrid,
                        test = aggregateTest)
      return(ci)
    }
  }

  # Oopps --------
  stop("Unsupported CI type!")
}

#' Retrieve Selection Adjusted P-Values
#'
#' @description A function for retrieving selection adjusted p-values
#' obtained after aggregate testing using the \code{\link{mvnQuadratic}},
#' \code{\link{mvnLinear}}, \code{\link{psatGLM}} or \code{\link{aggregatePvalues}} functions.
#'
#' @param object an object of class \code{\link{mvnQuadratic}},
#' \code{\link{glmQuadratic}}, or \code{\link{aggregatePvalues}}.
#'
#' @param type the type of p-value to compute, see the
#' the documentation of \code{\link{mvnQuadratic}} for details. If left
#' \code{NULL} then the first listed value in the \code{pvalue_type} field in
#' the original function call will be returned.
#'
#' @return a vector of p-values.
#'
getPval <- function(object, type = NULL) {
  # If aggregate pvalues --------------
  if(class(object)[1] == "aggregatePvalues") {
    return(object$p2C)
  } else if("mvnQuadratic" %in% class(object)) {
    return(getPvalQuadratic(object, type = type))
  } else if("mvnLinear" %in% class(object)) {
    return(getPvalLinear(object, type = type))
  } else {
    stop("Incorrect object type")
  }
}

getPvalQuadratic <- function(object, type = NULL) {
  if("psatGLM" %in% class(object)) {
    object$muhat <- object$betahat[-1]
    object$naiveMu <- object$naiveBeta[-1]
  }
  
  aggregateTest <- object$testType
  
  type <- type[1]
  if(is.null(type)) {
    type <- object$pvalue_type[1]
  }
  
  # Polyhedral ----------------------
  if(type == "polyhedral") {
    if(is.null(object$polyPval)) {
      return(getPolyCI(object$naiveMu, object$sigma, object$testMat,
                       object$threshold, confidence_level, FALSE,
                       test = aggregateTest, truncPmethod = truncPmethod)$pval)
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
                            object$threshold, confidence_level, FALSE,
                            test = aggregateTest)$pval
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

getPvalLinear <- function(object, type = NULL) {
  if("psatGLM" %in% class(object)) {
    object$muhat <- object$betahat[-1]
    object$naiveMu <- object$naiveBeta[-1]
  }
  
  aggregateTest <- object$testType
  
  type <- type[1]
  if(is.null(type)) {
    type <- object$pvalue_type[1]
  }
  
  # Polyhedral ----------------------
  if(type == "polyhedral") {
    if(is.null(object$polyPval)) {
      return(getPolyCI(object$naiveMu, object$sigma, object$contrast,
                       object$threshold, confidence_level, FALSE,
                       test = "linear")$pval)
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
                            object$threshold, confidence_level, FALSE,
                            test = aggregateTest)$pval
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

#' Retrieve Parameter Estimates from an mvnQuadratic Fit
#'
#' @description Retrieve parameter estimates from an
#' \code{\link{mvnQuadratic}} fit.
#'
#' @param object a \code{\link{mvnQuadratic}} object.
#'
#' @param type the type of estimates to return, should be
#' set to to either "mle" or "naive". If left \code{NULL} then
#' the mle will be returned.
#'
#' @return An estimate of the mean parameter of the normal
#' distribution.
coef.mvnQuadratic <- function(object, type = NULL, ...) {
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


#' Retrieve Parameter Estimates from an mvnLinear Fit
#'
#' @description Retrieve parameter estimates from an
#' \code{\link{mvnLinear}} fit.
#'
#' @param object a \code{\link{mvnLinear}} object.
#'
#' @param type the type of estimates to return, should be
#' set to to either "mle" or "naive". If left \code{NULL} then
#' the mle will be returned.
#'
#' @return An estimate of the mean parameter of the normal
#' distribution.
coef.mvnLinear <- function(object, type = NULL, ...) {
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


#' Retrieve Coefficient Estimates from a psatGLM Fit
#'
#' @description Retrieve coefficients estimates from a
#' \code{\link{psatGLM}} fit.
#'
#' @param object a \code{\link{psatGLM}} object.
#'
#' @param type the type of estimates to return, should be
#' set to to either "mle" or "naive". If left \code{NULL} then
#' the "mle" will be returned.
#'
#' @return An estimate of the regression coefficients of the
#' generalized linear model.
coef.psatGLM <- function(object, type = NULL, ...) {
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

predict.mvnQuadratic <- function(object, ...) {
  return(object$muhat)
}

#' Retrieve Linear Predictors from a psatGLM Fit
#'
#' @description Retrieve the linear predictor from a psatGLM
#' fit. This is the same as the \code{type = "link"} option in
#' \code{\link[stats]{predict.glm}}.
#'
#' @param object a \code{\link{psatGLM}} object.
#'
#' @param newX a new matrix of covariates. If left \code{NULL} the
#'  linear predictors for the data used for fitting the model will be
#'  returned.
#'
#' @return A vector of linear predictors.
#'
#' @seealso \code{\link{psatGLM}}
predict.psatGLM <- function(object, newX = NULL, ...) {
  if(is.null(newX)) {
    X <- object$X
  } else {
    X <- newX
  }

  X <- cbind(1, X)
  return(as.numeric(X %*% object$betahat))
}

#' Summarizing psatGLM Fits
#'
#' @description A summary method for \code{\link{psatGLM}} objects.
#'
#' @param object an object of type \code{\link{psatGLM}}.
#'
#' @param estimate_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param pvalue_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param ci_type see \code{\link{mvnQuadratic}} for details.
#'
#' @param confidence_level see \code{\link{mvnQuadratic}} for details.
#'
#' @details This is a summary method for summarizing the results of
#' post aggregate testing analysis with the \code{\link{psatGLM}}
#' method. The main output of the function is a table of regression coefficients
#' with confidence intervals and p-values computed using the desired methods.
#' This function is accompanied by a convenient printing function.
#'
#' @return A list containing a table of regression coefficients and details
#' regarding the inference methods used.
#'
#' @seealso \code{\link{psatGLM}}
summary.psatGLM <- function(object, estimate_type = NULL, pvalue_type = NULL,
                                 ci_type = NULL, confidence_level = NULL, ...) {
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
  sum$testType <- object$testType
  class(sum) <- "psatGLM_summary"
  return(sum)
}

print.psatGLM_summary <- function(x, ...) {
  sum <- x
  cat("Results for Post-Aggregate Testing Analysis \n \n")
  cat("Aggregate Test Type: ", sum$testType, "\n")
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






