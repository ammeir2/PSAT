#include <Rcpp.h>
using namespace Rcpp;

const double PIPI = 3.141592653589793238463 ;

void printVec(NumericVector x) { 
  for(int i = 0 ; i < x.length() ; i++) {
    Rcpp::Rcout<<x[i]<<" " ; 
  }
  Rcpp::Rcout<<"\n" ; 
  return ; 
}

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
double quadraticRobinsMonroe(NumericVector eta, 
                             NumericVector addVec, NumericVector multVec,
                             bool upper, double alpha,
                             double observed, double initU,
                             double threshold,
                             NumericMatrix sqrtTestMat,
                             NumericMatrix invSqrtTestMat,
                             NumericMatrix sampPrecision,
                             double stepSize,
                             int burnin, int mhiters, int rbIters) {
  // Initializing
  int p  = sqrtTestMat.nrow() ;
  double u = initU ;
  int tempBurnin = burnin ;
  NumericVector mu = NumericVector(p, 0.0) ;
  mu = addVec + multVec * u ;
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
    // printVec(mu) ; 
    //Rcpp::Rcout<<observed<<" "<<sum(samp * eta)<<"\n" ;
    if(upper) {
      c = k * (u - observed) / (i + 30) ;
      // c = k /  (i + 30) ;
      if(sum(samp * eta) > observed) {
        u -= c * alpha ;
      } else {
        u += c * (1 - alpha) ;
      }
    } else {
      c = k * (observed - u) / (i + 30) ;
      // c = k / (i + 30) ;
      if(sum(samp * eta) < observed) {
        u += c * alpha ;
      } else {
        u -= c * (1 - alpha) ;
      }
    }

    // Rcpp::Rcout<<samp[var]<<"\n" ;
    mu = addVec + multVec * u ;
    //printVec(mu) ; 
  }
  
  return u ;
}

// [[Rcpp::export]]
double linearRobinsMonroe(bool upper, double alpha,
                          double observed, 
                          double addVec, double multVec,
                          double initU,
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
    m = addVec + multVec * u ;
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
      c = k * double(u - observed) / double(i + 30) ;
      if(samp > observed) {
        u -= c * alpha ;
      } else {
        u += c * (1 - alpha) ;
      }
    } else {
      c = k * double(observed - u) / double(i + 30) ;
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

