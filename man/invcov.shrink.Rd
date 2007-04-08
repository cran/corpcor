\name{invcov.shrink}
\alias{invcov.shrink}
\alias{invcor.shrink}


\title{Fast Computation of the Inverse of the Covariance and Correlation Matrix}
\description{
  The functions \code{invcov.shrink} and \code{invcor.shrink} implement an
  algorithm to \emph{efficiently} compute 
  the inverses of shrinkage estimates of covariance (\code{\link{cov.shrink}}) 
  and correlation (\code{\link{cor.shrink}}).
}
\usage{
invcov.shrink(x, lambda, lambda.var, w, protect=0, verbose=TRUE)
invcor.shrink(x, lambda, w, protect=0, verbose=TRUE)
}
\arguments{
  \item{x}{a data matrix}
  \item{lambda}{the correlation shrinkage intensity (range 0-1).
                If \code{lambda} is not specified (the default) it is estimated
		using an analytic formula from Schaefer and Strimmer (2005) 
		- see \code{\link{cor.shrink}}.  
		For \code{lambda=0} the empirical correlations are recovered.}
  \item{lambda.var}{the variance shrinkage intensity (range 0-1). 
                If \code{lambda.var} is not specified (the default) it is estimated
	        using an analytic formula from Schaefer and Strimmer (2005) 
                - see \code{\link{var.shrink}}.  
	        For \code{lambda.var=0} the empirical variances are recovered.}	 
  \item{w}{optional: weights for each data point - if not specified uniform weights are assumed
            (\code{w = rep(1/n, n)} with \code{n = nrow(x)}).}	 
  \item{protect}{the fraction of correlation components protected against excessive individual 
                 component risk (default: 0, no protection by limiting translation)}	     
  \item{verbose}{output status while computing (default: TRUE)}
}
\details{
  The trick that allows the fast computation of the inverses  
   of the shrinkage covariance and correlation matrices is the
  Woodbury matrix identity - see, e.g., \url{http://en.wikipedia.org/wiki/Woodbury_matrix_identity}. 
The key insight from this identity  is that for inverting the covariance/correlation shrinkage estimator 
obtained from a n x p matrix you only need to invert a matrix of
the size of the rank of the data matrix (which in "small n, large p" setting
may mean substantial savings in computations). 

Note: at the moment the algorithmic savings are only available if \code{protect=0}.
}
\value{
  \code{invcov.shrink} returns the inverse of the output from \code{\link{cov.shrink}}. 
   
  \code{invcor.shrink} returns the inverse of the output from \code{\link{cor.shrink}}. 
}
\author{
  Juliane Schaefer (\url{http://www.stat.math.ethz.ch/~schaefer/}) 
  and Korbinian Strimmer (\url{http://strimmerlab.org}).
  
}
\references{

  Schaefer, J., and K. Strimmer. 2005.  A shrinkage approach to large-scale
  covariance estimation and implications for functional genomics. 
  Statist. Appl. Genet. Mol. Biol. \bold{4}:32.
  (\url{http://www.bepress.com/sagmb/vol4/iss1/art32/})
   
}
\seealso{\code{\link{cov.shrink}}, \code{\link{pcor.shrink}}, \code{\link{cor2pcor}}}

\examples{
\dontrun{

# load corpcor library
library("corpcor")

# generate data matrix
p <- 2000
n <- 10
X <- matrix(rnorm(n*p), nrow = n, ncol = p)

lambda <- 0.23  # some arbitrary lambda

# slow
system.time(
  W1 <-  solve(cov.shrink(X, lambda)) 
)

# very fast
system.time(
  W2 <- invcov.shrink(X, lambda)
)

# no difference
sum((W1-W2)^2)

}
}
\keyword{multivariate}
