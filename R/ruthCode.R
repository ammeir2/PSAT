#' Post Selection with Aggregate Testing for Independent P-Values
#'
#' @description Performs selection adjustment for independent p-values
#' following selection with an aggregate level test statistic.
#'
#' @param pmat An \code{MXn} matrix of original p-values , where \code{M} is
#' the number of features and \code{n} is the number of independent p-values
#' for each feature. If \code{globaltest="Pearson"}, the original p-values must
#' be one-sided p-values. No default.
#'
#' @param globaltest "Fisher" or "Pearson". The default is "Fisher". Set it to
#' "Pearson" only for two-sided alternatives, if common directionality is expected
#' for each feature.
#'
#' @param pthreshold The p-value selection threshold.
#' If \code{pthreshold=NULL}, only the global null p-value is computed. Must be
#' entered for computation of the conditional p-values. Either a scalar, or a vector
#' of length M with positive entries. See Details.
#'
#' @return List containing the global null p-value in \code{pF},  and the
#' conditional p-values in \code{p2C}. If a global null p-value is above
#' the threshold, then the conditional p-value vector corresponding to it is
#' a vector of NA's; if \code{pthreshold=NULL}, then the conditional p-value matrix
#' P2C is NA's.
#'
#' @details The Fisher test statisic is minus two times the sum of the log of
#' the input p-values, and the Fisher global null p-value is the probability
#' that a chi-squared distribution with \code{2n} degrees of freedom exceeds
#' the Fisher test statistc, where \code{n} is the number of input p-values.
#'
#' The Pearson test statistic is the maximum of the Fisher test statistic based
#' on left-sided p-values and the Fisher test statistic based on right-sided
#' p-values, and the Pearson global null p-value is twice the probability that
#' a chi-squared distribution with \code{2n} degrees of freedom exceeds the
#' Pearson test statistic. The global null p-value is output in pF. If the
#' hypotheses are one-sided, set globaltest= "Fisher"; if the hypotheses are
#' two-sided, set globaltest="Pearson". It is assumed that the p-values come from
#' continuous test statistics, so that left sided pvalue +  right sided pvalue = 1.
#' See Heller et al. (2016) for details about these global tests, about the
#' post-selection inference using these tests, and about extensions to other
#' global tests.
#'
#' The p-value threshold, pthreshold, is the threshold that the global null
#' p-value has to reach in order to be selected for post-selection inference.
#' For example,  for a bonferroni correction at level \code{alpha} for a family
#' of \code{m} global null hypotheses, it should be set to \code{alpha/m}. The
#' conditional p-values are computed only if \code{pF<=pthreshold}.
#'
#' @references Ruth Heller, Nilanjan Chatterjee, Abba Krieger, Jianxin Shi
#' (2016). Post-selection Inference Following Aggregate Level Hypothesis
#' Testing in Large Scale Genomic Data. bioRxiv:
#' \url{http://dx.doi.org/10.1101/058404}
#'
#' @example vignettes/aggregatePvalExamples.R
aggregatePvalues <- function(pmat, globaltest = c("Fisher", "Pearson"), pthreshold = NULL){

  if(sum(is.na(pmat)) > 0){
    warning("NAs in the conditional pvalue compuations")
  }
  if(!length(globaltest)) {
    stop("globaltest must be Fisher or Pearson")
  }
  globaltest <- globaltest[1]

  M = dim(pmat)[1]
  if(length(pthreshold) == 1) {
    pthreshold <- rep(pthreshold, M)
  }


  p2C <- matrix(NA, nrow = M, ncol = dim(pmat)[2])
  pF <- rep(NA, M)

  if(!(globaltest %in% c("Fisher", "Pearson"))) {
    stop("globaltest must be Fisher or Pearson")
  }

  for(i in 1:M) {
    p <- as.double(pmat[i, ])
    N <- length(p)
    whichNotNA <- which(!is.na(p))
    p <- p[!is.na(p)]
    n <- length(p)
    if(n > 0) {
      if(globaltest == "Fisher") {
        tp <- -2 * log(p)
        pF[i] <- 1 - pchisq(sum(tp), 2 * n)
        if(!is.null(pthreshold)) {
          if(length(pthreshold) < M) {
            stop("threshold t for post-selection computation does not match the number of rows in pmat")
          }
          t <- qchisq(1-pthreshold[i], 2*n)

          if(sum(tp) >= t){
            p2c <- rep(NA, n)
            iscond <- rep(0, n)
            p2cond <- prod(p) / exp(-t/2)
            if(prod(sort(p)[2:n]) <= exp(-t/2)){ #no payment for row selection
              p2C[i, ] <- p
            } else {
              for(j in 1:n) {
                iscond[j] <- sum(sum(tp[-j]) < t)
                p2c[j] = ifelse(iscond[j] == 1, p2cond, p[j])
              }#for
              p2C[i, whichNotNA] <- p2c
            }
          }#end    if ( sum(tp)>=t)
        }#end    if (!is.null(pthreshold))
      }# end if (globaltest=="Fisher")


      if(globaltest == "Pearson"){
        tpL <- -2 * log(p)
        tpR <- -2 * log(1-p)
        QL <- sum(tpL)
        QR <- sum(tpR)
        QC <- max(QL, QR)
        pF[i] <- 2 * (1 - pchisq(QC, 2 * n))

        if(!is.null(pthreshold)) {
          if(length(pthreshold) < M) {
            stop("threshold t for post-selection computation does not match the number of rows in pmat")
          }
          t <-  qchisq(1 - pthreshold[i] / 2, 2 * n)
          if(QC >= t){
            p2c <- rep(NA, n)
            numc <- rep(NA, n)
            denomc <- rep(NA, n)
            for(j in 1:n){
              logAL <- -t/2 - sum(log(p[-j]))
              AL <- exp(logAL)
              logAR <- -t/2 - sum(log((1-p)[-j]))
              AR <- exp(logAR)
              p2sided <- 2 * min(p[j], 1 - p[j])
              if(AL + AR >= 1){
                p2c[j] <- p2sided
              }
              if(AL+AR < 1) {
                denomc[j] <- AL + AR
                numc[j] <- min(p2sided / 2, AL) + min(p2sided / 2, AR) +
                  max(AL - (1 - p2sided / 2), 0) + max(AR - (1 - p2sided/2), 0)
                p2c[j] <- numc[j] / denomc[j]
              }
            }#for
            p2C[i, whichNotNA] <- p2c
          }#end if (QC>=t)
        }#end if (!is.null(pthreshold))
      }# end if (globaltest=="Pearson")
    }# end if (n>0)
  }#end for (i in 1:M)

  result <- list(pF = pF, p2C = p2C)
  class(result) <- "aggregatePvalues"
  return(result)
}#Conditional.pvalues
