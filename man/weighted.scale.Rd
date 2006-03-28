\name{weighted.scale}
\alias{weighted.scale}
\alias{weighted.var}
\alias{weighted.moments}

\title{Weighted Expectations and Variance}

\description{
  \code{weighted.var} estimate the unbiased variance taking into account 
  data weights.

  \code{weighted.moments} produces the weighted mean and weighted variance
  for each column of a matrix.
  
  \code{weighted.scale} centers and standardized a matrix using 
  the weighted mean and variance.
 }

\usage{
  weighted.var(xvec, w) 
  weighted.moments(x, w)
  weighted.scale(x, w, center=TRUE, scale=TRUE, wm)
}

\arguments{
  \item{xvec}{a vector}
  \item{x}{a matrix}
  \item{w}{data weights}
  \item{center}{logical value}
  \item{scale}{logical value}
  \item{wm}{list with variances and means as returned by \code{weighted.moments}.
            Will be computed from \code{x} if missing.}
}


\value{
  A rescaled matrix (\code{weighted.scale}), a list containg the column means and 
  variances (\code{weighted.moments}), or simple number (\code{weighted.var})
}


\author{
  Rainer Opgen-Rhein (\url{http://www.statistik.lmu.de/~opgen/}) and
  Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
}

\seealso{\code{\link{weighted.mean}}.}

\examples{
# load corpcor library
library("corpcor")

# generate some data
p <- 5
n <- 5
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
w <- c(1,1,1,3,3)/9


# standardize matrix
scale(X)
weighted.scale(X)
weighted.scale(X, w) # take into account data weights

}
\keyword{multivariate}
