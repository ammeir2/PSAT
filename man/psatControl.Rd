% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalInference.R
\name{psatControl}
\alias{psatControl}
\title{Creates a list of parameters for use with PSAT inference functions}
\usage{
psatControl(switchTune = NULL, nullMethod = c("RB", "zero-quantile"),
  nSamples = NULL, sgdStep = NULL, nsteps = 1000, trueHybrid = FALSE,
  rbIters = NULL, optimMethod = c("Nelder-Mead", "SGD"),
  truncPmethod = c("UMPU", "symmetric"), quadraticSampler = c("tmg",
  "PSAT"))
}
\arguments{
\item{switchTune}{tuning parameter for regime switching confidence intervals, 
should be a number between 0 and confidence_level.}

\item{nullMethod}{method for compute global-null confidence intervals. Robins-Monroe (RB) should be used.}

\item{nSamples}{number of samples to be taken from the null distribution.}

\item{sgdStep}{number of stochastic gradient steps to take, only applicable for 
non-wald quadratic tests when \code{optimMethod} "SGD" is used.}

\item{trueHybrid}{whether Robins-Monroe computation should be performed for computing 
confidence intervals for all confidence intervals or only when they may improve power.}

\item{rbIters}{number of steps to take when computing confidence intervals with the Robins-Monroe
procedure.}

\item{optimMethod}{optimization method to be used when computing the conditional MLE in
in inference after testing with a non-wald aggregate test. Nelder-Mead as implemented in the 
\code{\link[stats]{optim}} function.}

\item{truncPmethod}{the type of test to use when computing the polyhedral p-values,
options are either the UMPU test or a symmetric test.}

\item{quadraticSampler}{which quadratic sampler to use? Choices are either the 
Hamiltionian Montel-Carlo method implemented in \code{\link[tmg]{rtmg}},
or a Gibbs sampler implemented in this package.}
}
\description{
Creates a list with additional parameters for use with 
\code{\link{mvnQuadratic}}, \code{\link[psat]{mvnLinear}}, and \code{\link[PSAT]{psatGLM}}.
}
