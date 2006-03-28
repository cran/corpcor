\name{pcor.shrink}
\alias{pcor.shrink}
\alias{pcov.shrink}

\title{Shrinkage Estimates of Partial Correlation and Partial Covariance}
\description{
  The functions \code{pcor.shrink} and \code{pcov.shrink} provide shrinkage estimates
  the matrix of partial correlations and partical covariance, respectively.  
}
\usage{
pcor.shrink(x, lambda, w, verbose=TRUE)
pcov.shrink(x, lambda, lambda.var, w, verbose=TRUE)
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
  \item{w}{optional: weights for each data point - if not specified uniform weights
            are assumed (\code{w = rep(1/n, n)} with \code{n = dim(x)[1]}).}	 
  \item{verbose}{report progress while computing (default: TRUE)}
}
\details{
   The partial covariance is simply the inverse of the covariance matrix, with
   the sign of the off-diagonal elements reversed.  The partial correlation matrix
   is the standardized partial covariance matrix.
  
   Note that \code{pcor.shrink(x)} is numerically the same as \code{cor2pcor(cor.shrink(x))}.
   However, \code{pcor.shrink} is computationally *much* faster.
   
   For details about the shrinkage procedure (from Schaefer and Strimmer (2005))
   please consult the help page of \code{\link{cov.shrink}}.
   
}
\value{
  \code{pcor.shrink} returns the partical correlation matrix.
  
  \code{pcov.shrink} returns the partial covariance matrix. 
}
\author{
  Juliane Schaefer (\url{http://www.stat.math.ethz.ch/~schaefer/}) 
  and Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
  
}
\references{

  Schaefer, J., and Strimmer, K. (2005).  A shrinkage approach to large-scale
  covariance estimation and implications for functional genomics. 
  Statist. Appl. Genet. Mol. Biol.\bold{4}:32.
  (\url{http://www.bepress.com/sagmb/vol4/iss1/art32/})  
}
\seealso{\code{\link{cov.shrink}}, \code{\link{cor2pcor}}}

\examples{
# load corpcor library
library("corpcor")

# generate data matrix
p <- 50
n <- 10
X <- matrix(rnorm(n*p), nrow = n, ncol = p)


# partial covariance
pco <- pcov.shrink(X)

# partial correlations (fast and recommend way)
pcr1 <- pcor.shrink(X)

# other possibilites to estimate partial correlations
pcr2 <- cor2pcor( cor.shrink(X) )
pcr3 <- cov2cor( pco ) # standarize partial covariance


# all the same
sum((pcr1 - pcr2)^2)
sum((pcr2 - pcr3)^2)
sum((pcr3 - pcr1)^2)

}
\keyword{multivariate}
