\name{varcov}
\alias{varcov}

\title{Variance of the Entries of the Covariance Matrix}
\description{
  
  \code{varcov} computes both the standard empirical covariance matrix (\eqn{S}) and  a matrix 
   containing the estimated variances \eqn{ var(s_{ij})} of the individual components of (\eqn{S}).
 
}
\usage{
varcov(x, type=c("unbiased", "ML"), verbose=TRUE)}
\arguments{
  \item{x}{a data matrix}
  \item{type}{estimate unbiased or ML variant of empirical covariance matrix
              (default: "unbiased", as in \code{\link{cov}})}
  \item{verbose}{report progress while computing (default: TRUE)}
	      
}
\details{
   The variances of the empirical covariance matrix are needed to
   compute shrinkage estimates - see Schaefer and Strimmer (2005). 
}

\value{ 
  \code{varcov} returns a list with the following two components:
  \item{S}{ the empirical covariance matrix.   } 
  \item{var.S}{ a matrix containing the variances of each element in S.}
}

\author{
  Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
  
}

\references{

  Schaefer, J., and Strimmer, K. (2005).  A shrinkage approach to large-scale
  covariance estimation and implications for functional genomics.  Submitted
  to SAGMB.  
}

\seealso{\code{\link{cov.shrink}}.}

\examples{
# load corpcor library
library("corpcor")

# small n, large p
p <- 100
n <- 20

# generate random pxp covariance matrix
sigma <- matrix(rnorm(p*p),ncol=p)
sigma <- crossprod(sigma)+ diag(rep(0.1, p))

# simulate multinormal data of sample size n  
sigsvd <- svd(sigma)
Y <- t(sigsvd$v \%*\% (t(sigsvd$u) * sqrt(sigsvd$d)))
X <- matrix(rnorm(n * ncol(sigma)), nrow = n) \%*\% Y


# estimate covariance matrix
S <- cov(X)

# varcov produces the same results as cov
vc <- varcov(X)
sum(abs(vc$S-S))
}
\keyword{multivariate}
