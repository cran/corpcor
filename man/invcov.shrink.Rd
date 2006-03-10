\name{invcov.shrink}
\alias{invcov.shrink}
\alias{invcor.shrink}


\title{Fast Computation of the Inverse of the Covariance and Correlation Matrix}
\description{
  The functions \code{invcov.shrink} and \code{invcor.shrink} implement an
  algrithm to *efficiently* compute 
  the inverses of shrinkage estimates of covariance (\code{\link{cov.shrink}}) 
  and correlation (\code{\link{cor.shrink}}).
}
\usage{
invcov.shrink(x, lambda, w, verbose=TRUE)
invcor.shrink(x, lambda, w, verbose=TRUE)
}
\arguments{
  \item{x}{a data matrix}
  \item{lambda}{the shrinkage intensity (range 0-1). If \eqn{\lambda} is is not specified (the default)
        a suitable value is automatically chosen such that the resulting shrinkage estimate
	 has minimal MSE (see below for details).  For \eqn{\lambda=0} not shrinkage occurs
	 and the empirical estimates are computed.}
  \item{w}{optional: weights for each data point - if not specified uniform weights are assumed
            (\code{w = rep(1/n, n)} with \code{n = dim(x)[1]}).}	 
  \item{verbose}{output status while computing (default: TRUE)}
}
\details{
  The trick that allows the fast computation of the inverses  
   of the shrinkage covariance and correlation matrices is the
  Woodbury matrix identity (\url{http://en.wikipedia.org/wiki/Woodbury_matrix_identity}).
}
\value{
  \code{invcov.shrink} returns the inverse of the output from \code{\link{cov.shrink}}. 
   
  \code{invcor.shrink} returns the inverse of the output from \code{\link{cor.shrink}}. 
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
