\name{cov.bagged}
\alias{cov.bagged}
\alias{cor.bagged}
\alias{pcor.bagged}

\title{Bagged Estimates of Covariance and (Partial) Correlation}

\usage{
cov.bagged(x, R=1000, \dots)
cor.bagged(x, R=1000, \dots)
pcor.bagged(x, R=1000, \dots)
}

\arguments{
  \item{x}{data matrix or data frame}
  \item{R}{number of bootstrap replicates (default: 1000)}
  \item{\dots}{options passed to \code{\link{cov}}, 
               \code{\link{cor}}, and \code{\link{pseudoinverse}}  }
}
\description{
  \code{cov.bagged}, \code{cor.bagged}, and \code{pcor.bagged} calculate
  the bootstrap aggregated (=bagged) versions of the covariance and
  (partial) covariance estimators.  
}

\details{
  Bootstrap aggregation, or ``bagging'', was first suggested by Breiman (1996)
  as a means to improve an estimator using the bootstrap. The bagged estimate 
  corresponds simply to the mean of the bootstrap sampling distribution.
  
  Bagging is essentially a non-parametric variance reduction method.
  In Schaefer and Strimmer (2005a,b) the inverse of the bagged correlation matrix
  is used to estimate graphical Gaussian models from sparse microarray data.

  Note that bagging is computatationally quite expensive.  
  See  \code{\link{cov.shrink}} for alternative (and better) approach 
  to variance-reduced covariance estimation.  
}
\value{
  A symmetric matrix.
}
\author{
  Juliane Schaefer (\url{http://www.stat.math.ethz.ch/~schaefer/}) and
  Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
}
\references{
  Breiman, L. (1996). Bagging predictors. \emph{Machine Learning}, \bold{24}, 123--140.

  Schaefer, J., and Strimmer, K. (2005a).  An empirical Bayes approach to inferring
  large-scale gene association networks. \emph{Bioinformatics} \bold{21}:754-764.
  
  Schaefer, J., and Strimmer, K. (2005b). Learning large-scale graphical Gaussian 
  models from genomic data. Proceedings of CNET 2004, Aveiro, Pt. (AIP)
  
  
}
\seealso{\code{\link{cov.shrink}}.}

\examples{
# load corpcor library
library("corpcor")

# some statistics on the US states
data(state)
us.states <- t(state.x77)
dim(us.states) # sample size: 8, number of variables: 50

# estimates of correlation
c <- cor(us.states)
bc <- cor.bagged(us.states)
sc <- cor.shrink(us.states)

# positive definiteness of bagged correlation and shrinkage correlation
is.positive.definite(c)
is.positive.definite(bc)
is.positive.definite(sc)

# rank and condition
rank.condition(c)
rank.condition(bc)
rank.condition(sc) # overall best!

}
\keyword{multivariate}
