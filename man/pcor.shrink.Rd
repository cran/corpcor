\name{pcor.shrink}
\alias{pcor.shrink}
\alias{pvar.shrink}

\title{Shrinkage Estimates of Partial Correlation and Partial Variance}
\description{
  The functions \code{pcor.shrink} and \code{pvar.shrink} compute shrinkage estimates
  of partial correlation and partial variance, respectively.  
}
\usage{
pcor.shrink(x, lambda, w, protect=0.01, verbose=TRUE)
pvar.shrink(x, lambda, lambda.var, w, protect=0.01, verbose=TRUE)
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
            are assumed (\code{w = rep(1/n, n)} with \code{n = nrow(x)}).}	 
  \item{protect}{the fraction of correlation components protected against excessive individual 
                 component risk (default: 0.01)}	     
  \item{verbose}{report progress while computing (default: TRUE)}
}
\details{
   The partial variance \eqn{var(X_k | rest)} is the variance of \eqn{X_k} conditioned on the
   remaining variables.  It equals the inverse of the corresponding diagonal entry
   of the precision matrix (see Whittaker 1990).
   
   The partial correlations \eqn{corr(X_k, X_l | rest)} is the correlation between
   \eqn{X_k} and \eqn{X_l} conditioned on the remaining variables. It equals the sign-reversed
   entries of the off-diagonal entries of the precision matrix, standardized
   by the the squared root of the associated inverse partial variances.
  
  
   Note that using \code{pcor.shrink(x)} \emph{much} faster than 
   \code{cor2pcor(cor.shrink(x))}.
   
   For details about the shrinkage procedure consult Schaefer and Strimmer (2005)
   and the help page of \code{\link{cov.shrink}}.
   
}
\value{
  \code{pcor.shrink} returns the partial correlation matrix.
  
  \code{pvar.shrink} returns the partial variances. 
}
\author{
  Juliane Schaefer (\url{http://www.stat.math.ethz.ch/~schaefer/}) 
  and Korbinian Strimmer (\url{http://strimmerlab.org}).
  
}
\references{

  Schaefer, J., and Strimmer, K. (2005).  A shrinkage approach to large-scale
  covariance estimation and implications for functional genomics. 
  Statist. Appl. Genet. Mol. Biol.\bold{4}:32.
  (\url{http://www.bepress.com/sagmb/vol4/iss1/art32/})  
  
  Whittaker J. (1990).  Graphical Models in Applied Multivariate Statistics.
  John Wiley, Chichester.
    
}
\seealso{\code{\link{invcov.shrink}}, \code{\link{cov.shrink}}, \code{\link{cor2pcor}}}

\examples{
# load corpcor library
library("corpcor")

# generate data matrix
p <- 50
n <- 10
X <- matrix(rnorm(n*p), nrow = n, ncol = p)


# partial variance
pv <- pvar.shrink(X)
pv

# partial correlations (fast and recommend way)
pcr1 <- pcor.shrink(X)

# other possibilities to estimate partial correlations
pcr2 <- cor2pcor( cor.shrink(X) )


# all the same
sum((pcr1 - pcr2)^2)

}
\keyword{multivariate}
