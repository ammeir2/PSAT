aggregatePvalues <- function(pmat, globaltest = c("Fisher", "Pearson"), pv.threshold = NULL){

  if(sum(is.na(pmat)) > 0){
    warning("NAs in the conditional pvalue compuations")
  }
  if(!length(globaltest)) {
    stop("globaltest must be Fisher or Pearson")
  }
  globaltest <- globaltest[1]

  M = dim(pmat)[1]
  if(length(pv.threshold) == 1) {
    pv.threshold <- rep(pv.threshold, M)
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
        if(!is.null(pv.threshold)) {
          if(length(pv.threshold) < M) {
            stop("threshold t for post-selection computation does not match the number of rows in pmat")
          }
          t <- qchisq(1-pv.threshold[i], 2*n)

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
        }#end    if (!is.null(pv.threshold))
      }# end if (globaltest=="Fisher")


      if(globaltest == "Pearson"){
        tpL <- -2 * log(p)
        tpR <- -2 * log(1-p)
        QL <- sum(tpL)
        QR <- sum(tpR)
        QC <- max(QL, QR)
        pF[i] <- 2 * (1 - pchisq(QC, 2 * n))

        if(!is.null(pv.threshold)) {
          if(length(pv.threshold) < M) {
            stop("threshold t for post-selection computation does not match the number of rows in pmat")
          }
          t <-  qchisq(1 - pv.threshold[i] / 2, 2 * n)
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
        }#end if (!is.null(pv.threshold))
      }# end if (globaltest=="Pearson")
    }# end if (n>0)
  }#end for (i in 1:M)

  result <- list(pF = pF, p2C = p2C)
  class(result) <- "aggregatePvalues"
  return(result)
}#Conditional.pvalues
