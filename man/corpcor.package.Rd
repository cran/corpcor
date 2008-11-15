\encoding{latin1}
\name{corpcor-package}
\alias{corpcor-package}
\docType{package}
\title{The corpcor package}

\description{
This package implements an analytic shrinkage approach
  for inferring the covariance matrix.  The estimator is statistically
  highly accurate and efficient, applicable to "small n, large p" data, 
  and always returns a positive definite and well-conditioned matrix. 
  Nevertheless, this method requires only little a priori modeling and is
  computationally cheap. In addition to covariance estimation the
  package contains similar functions for inferring variances, correlations,
  partial correlations, partial variances, and regression coefficients.  
  Furthermore, it provides  functions for fast SVD computation, for computing 
  the pseudoinverse, for checking the rank and positive definiteness of a
  matrix, and for the computationally fast inversion of the covariance 
  and correlation matrix.
}

\author{Juliane Sch\"afer, Rainer Opgen-Rhein, and Korbinian Strimmer (\url{http://strimmerlab.org/})}
\references{
  See website: \url{http://strimmerlab.org/software/corpcor/}
}
\keyword{multivariate}
\seealso{
\code{\link{cov.shrink}, \link{invcov.shrink}, \link{pcor.shrink}, 
 \link{mvr.shrink}, \link{fast.svd}.}
}
