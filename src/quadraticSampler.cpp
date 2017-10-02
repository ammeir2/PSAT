#include <Rcpp.h>
using namespace Rcpp;

const double PIPI = 3.141592653589793238463 ;

double sampleUnivTruncNorm(double mu, double sd, double threshold) {
  double u = runif(1)[0] ;
  double phiThreshold, sample ;

  if(threshold > 0) {
    phiThreshold = R::pnorm5(threshold, mu, sd, 0, 0) ;
    sample = R::qnorm5(u * phiThreshold, mu, sd, 0, 0) ;
  } else {
    phiThreshold = R::pnorm5(threshold, mu, sd, 1, 0) ;
    sample = R::qnorm5(u * phiThreshold, mu, sd, 1, 0) ;
  }

  return sample ;
}

double sampleUnivTruncNorm(double mu, double sd, double threshold, bool positive) {
  double u = runif(1)[0] ;
  double phiThreshold, sample ;

  if(positive) {
    phiThreshold = R::pnorm5(threshold, mu, sd, 0, 0) ;
    sample = R::qnorm5(u * phiThreshold, mu, sd, 0, 0) ;
  } else {
    phiThreshold = R::pnorm5(threshold, mu, sd, 1, 0) ;
    sample = R::qnorm5(u * phiThreshold, mu, sd, 1, 0) ;
  }

  return sample ;
}


double computeConditionalMean(NumericVector mu, NumericVector samp,
                              NumericMatrix precision, int index) {
  double result = 0 ;

  for(int j = 0; j < mu.length() ; j ++) {
    if(j != index) {
      result += precision(index, j) * (samp[j] - mu[j]) ;
    }
  }

  result = result / precision(index, index) ;
  result = mu[index] - result ;
  return result ;
}

double computeSquaredNorm(NumericVector x, int index) {
  double result = 0 ;

  for(int i = 0 ; i < x.length() ; i++) {
    if(i != index) {
      result += std::pow(x[i], 2.0) ;
    }
  }

  return result ;
}

void voidSampleQuadraticMat(NumericMatrix sampleMatrix,
                            NumericVector init, NumericVector mean,
                            NumericMatrix precision, double threshold,
                            int burnin, int trim) {
  int nSamples = sampleMatrix.nrow() ;
  int p = mean.length() ;
  int sampNum = 0;
  NumericVector sample = clone(init) ;
  double sd, condMean ;
  double normMinus, currentThreshold;
  double pplus, pminus, unif ;
  burnin = std::max(burnin, trim) ;
  int iteration = 0 ;

  while(sampNum < nSamples) {
    iteration++ ;
    for(int j = 0; j < p ; j++) {
      sd = 1 / sqrt(precision(j, j)) ;
      condMean = computeConditionalMean(mean, sample, precision, j) ;
      normMinus = computeSquaredNorm(sample, j) ;
      currentThreshold = threshold - normMinus ;
      if(currentThreshold <= 0) {
        sample[j] = rnorm(1, condMean, sd)[0] ;
      } else {
        currentThreshold = std::sqrt(currentThreshold) ;
        pplus = R::pnorm5(currentThreshold, condMean, sd, 0, 1) ;
        pminus = R::pnorm5(-currentThreshold, condMean, sd, 1, 1) ;
        pplus = 1 / (1 + std::exp(pminus - pplus)) ;
        unif = runif(1)[0] ;
        if(unif < pplus) {
          sample[j] = sampleUnivTruncNorm(condMean, sd, currentThreshold) ;
        } else {
          sample[j] = sampleUnivTruncNorm(condMean, sd, -currentThreshold) ;
        }
      }
    }

    if(iteration >= burnin && (iteration % trim == 0)) {
      sampleMatrix(sampNum++, _) = sample ;
    }
  }
}

// [[Rcpp::export]]
NumericMatrix sampleQuadraticMVTcpp(NumericVector init, NumericVector mean,
                                    NumericMatrix precision, double threshold,
                                    int nSamples, int burnin, int trim) {
  NumericMatrix sampleMatrix(nSamples, init.length()) ;
  voidSampleQuadraticMat(sampleMatrix, init,  mean,
                         precision, threshold, burnin, trim) ;
  return sampleMatrix ;
}

void innerProduct(NumericVector result, NumericMatrix x, NumericVector y) {
  int n = x.nrow() ;
  int p = x.ncol() ;

  for(int i = 0; i < n ; i++) {
    result[i] = 0 ;
    for(int j = 0; j < p ; j++) {
      result[i] += x(i, j) * y[j] ;
    }
  }

}

// [[Rcpp::export]]
double quadraticRobinsMonroe(int var, bool upper, double alpha,
                                    double observed, double initU,
                                    double threshold,
                                    NumericMatrix sqrtTestMat,
                                    NumericMatrix invSqrtTestMat,
                                    NumericMatrix sampPrecision,
                                    double stepSize,
                                    int burnin, int mhiters, int rbIters) {
  // Initializing
  var -= 1;
  int p  = sqrtTestMat.nrow() ;
  double u = initU ;
  int tempBurnin = burnin ;
  NumericVector mu = NumericVector(p, 0.0) ;
  mu[var] = u ;
  NumericVector sampMu = clone(mu) ;
  NumericVector samp = clone(mu) ;
  NumericMatrix rawsamp = NumericMatrix(1, p) ;
  innerProduct(rawsamp(0, _), sqrtTestMat, mu) ;
  // Step Size ;
  // double qz = R::qnorm5(1 - 0.01, 0, 1, 1, 0) ;
  double k = stepSize ;
  double c ;

  for(int i = 1 ; i < rbIters + 1 ; i++) {
    if(i == 1) {
      tempBurnin = mhiters ;
    }

    // Sampling
    innerProduct(sampMu, sqrtTestMat, mu) ;
    voidSampleQuadraticMat(rawsamp, rawsamp(0, _),
                           sampMu, sampPrecision, threshold,
                           tempBurnin, 1) ;
    innerProduct(samp, invSqrtTestMat, rawsamp(0, _)) ;

    // Taking step
    if(upper) {
      c = k * (u - observed) / (i + 30) ;
      // c = k /  (i + 30) ;
      if(samp[var] > observed) {
        // Rcpp::Rcout<<k<<" "<<samp[var]<<" " ;
        u -= c * alpha ;
      } else {
        u += c * (1 - alpha) ;
      }
    } else {
      c = k * (observed - u) / (i + 30) ;
      // c = k / (i + 30) ;
      if(samp[var] < observed) {
        u += c * alpha ;
      } else {
        u -= c * (1 - alpha) ;
      }
    }

    // Rcpp::Rcout<<samp[var]<<"\n" ;
    mu[var] = u ;
  }

  return u ;
}

// [[Rcpp::export]]
double linearRobinsMonroe(bool upper, double alpha,
                          double observed, double initU,
                          double contrastCoordinate,
                          NumericVector threshold,
                          double contsd,
                          double condSD,
                          double regConst,
                          double stepSize,
                          int rbIters) {
  // Initializing
  double u = initU ;
  double k = stepSize ;
  double c, m, condMu;
  double posProb, negProb ;
  double contSamp ;
  double samp ;

  for(int i = 1 ; i < rbIters + 1 ; i++) {
    // Sampling contrast
    m = contrastCoordinate * u;
    negProb = R::pnorm5(threshold[1], m, contsd, 1, 1) ;
    posProb = R::pnorm5(threshold[2], m, contsd, 0, 1) ;
    negProb = 1 / (1 + std::exp(posProb - negProb)) ;
    if(runif(1)[0] < negProb) {
      contSamp = sampleUnivTruncNorm(m, contsd, threshold[0], false) ;
    } else {
      contSamp = sampleUnivTruncNorm(m, contsd, threshold[1], true) ;
    }
    // Rcpp::Rcout<<contSamp<<" " ;

    // Sampling Coordinate
    condMu = u + regConst * (contSamp - m) ;
    samp = R::rnorm(condMu, condSD) ;

    // Taking step
    if(upper) {
      c = k * (u - observed) / (i + 30) ;
      if(samp > observed) {
        u -= c * alpha ;
      } else {
        u += c * (1 - alpha) ;
      }
    } else {
      c = k * (observed - u) / (i + 30) ;
      if(samp < observed) {
        u += c * alpha ;
      } else {
        u -= c * (1 - alpha) ;
      }
    }
  }

  // Rcpp::Rcout<<"\n";
  return u ;
}


double innerProduct(NumericVector x, NumericVector y) {
  double result = 0 ;
  for(int i = 0; i < x.length() ; i++) {
    result += x[i] * y[i] ;
  }

  return result;
}

void vectorAddition(NumericVector result, NumericVector x,
                    NumericVector y, double multiplier) {
  double p = result.length() ;
  for(int j = 0; j < p ; j++) {
    result[j] = x[j] + multiplier * y[j] ;
  }
}

// NumericMatrix sampleLogisticCpp(NumericVector y,
//                                 NumericMatrix X,
//                                 NumericVector Xy,
//                                 NumericVector Xp,
//                                 NumericVector probs,
//                                 NumericVector beta,
//                                 double intercept,
//                                 NumericMatrix testMat,
//                                 double threshold,
//                                 int nSamples, int burnin, int trim,
//                                 NumericVector returnY) {
//   // Declarations and preliminaries
//   int n = y.length() ;
//   int p = beta.length() ;
//   double meanProb ;
//   // NumericMatrix sampleMatrix(nSamples, p + 1) ;
//   NumericMatrix sampleMatrix(nSamples, y.length()) ;
//   NumericVector newXy = clone(Xy) ;
//   NumericVector weights = probs * (1 - probs) ;
//   NumericVector innerP(p) ;
//   NumericVector diff(p) ;
//   int samp ;
//   double newy, MHratio, testStat;
//   double prob, sumy, newSumy ;
//   sumy = sum(y) ;
//
//   // Making sure changes are not too unlikely
//   double minProb = min(probs) ;
//   for(int k = 0 ; k < 2 ; k ++) {
//     for(int i = 0; i < n ; i ++) {
//       if(k == 0) {
//         probs[i] -= minProb ;
//         probs[i] = std::max(probs[i], 0.05) ;
//       } else {
//         probs[i] = std::min(probs[i] / meanProb * 0.6, 0.95) ;
//       }
//       meanProb = mean(probs) ;
//     }
//   }
//
//   // Sampling!
//   int iter = 0 ;
//   int m = 0 ;
//   while(m < sampleMatrix.nrow()) {
//     iter ++ ;
//     for(int i = 0; i < n ; i++) {
//       samp = i ;
//       newy = R::rbinom(1, probs[samp]);
//       // If newy == oldy then no change is necessary
//       if(newy == y[samp]) continue ;
//
//       // Creating Xy proposal
//       if(newy == 0) {
//         vectorAddition(newXy, Xy, X(samp, _), -1) ;
//         newSumy = sumy - 1 ;
//         prob = 1 - probs[samp] ;
//       } else {
//         vectorAddition(newXy, Xy, X(samp, _), 1) ;
//         newSumy = sumy + 1 ;
//         prob = probs[samp] ;
//       }
//
//       // Do MHratio
//       MHratio = innerProduct(newXy, beta) - innerProduct(Xy, beta) ;
//       MHratio += intercept * (newSumy - sumy) ;
//       MHratio += -log(prob) + log(1 - prob) ;
//       MHratio = std::exp(MHratio) ;
//       //Rcpp::Rcout<<newSumy<<" ";
//       if(runif(1)[0] < MHratio) {
//         // Check that constraint is satisfied
//         vectorAddition(diff, newXy, Xp, -1) ;
//         innerProduct(innerP, testMat, diff) ;
//         testStat = innerProduct(diff, innerP) ;
//         if(testStat > threshold) {
//           //Rcpp::Rcout<<testStat<<"\n" ;
//           y[samp] = newy ;
//           sumy = newSumy ;
//           vectorAddition(Xy, newXy, newXy, 0) ;
//           continue ;
//         }
//       }
//
//       // If proposal isn't good, revert to previous
//       vectorAddition(newXy, Xy, Xy, 0) ;
//     }
//     //Rcpp::Rcout<<"\n ";
//
//     if((iter > burnin) & (iter % trim == 0)) {
//       // for(int j = 0 ; j < p ; j++) {
//       //   sampleMatrix(m, j + 1) = Xy[j] ;
//       // }
//       // sampleMatrix(m++, 0) = sumy ;
//       sampleMatrix(m++, _) = y ;
//     }
//   }
//
//   returnY = y ;
//   return sampleMatrix ;
// }
//
