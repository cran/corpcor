\name{cor.shrink}
\alias{cov.shrink}
\alias{cor.shrink}
\alias{pcor.shrink}

\title{Shrinkage Estimates of Covariance and (Partial) Correlation}
\description{
  The functions \code{cov.shrink}, \code{cor.shrink}, and \code{pcor.shrink} implement a shrinkage approach 
  to estimate covariance and (partial) correlation matrices, as suggested in Schaefer and Strimmer (2005).
  
  The advantages of using this approach in comparison with the standard empirical estimates 
  (\code{\link{cov}} and \code{\link{cor}}) are that the shrinkage estimates 
  \enumerate{
    \item are always  positive definite,
    \item  well conditioned (so that the inverse always exists), and 
    \item exhibit (sometimes dramatically) better mean squared error.}
  Furthermore, they are inexpensive to compute and do not require any
  tuning parameters (the shrinkage intensity is analytically estimated from the data).
  
}
\usage{
cov.shrink(x, lambda, verbose=TRUE)
cor.shrink(x, lambda, verbose=TRUE)
pcor.shrink(x, lambda, verbose=TRUE)
}
\arguments{
  \item{x}{a data matrix}
  \item{lambda}{the shrinkage intensity (range 0-1). If \eqn{\lambda} is is not specified (the default)
        a suitable value is automatically chosen such that the resulting shrinkage estimate
	 has minimal MSE (see below for details).}
  \item{verbose}{report progress while computing (default: TRUE)}
}
\details{
   
  \code{cor.shrink} computes a shrinkage estimate \eqn{R^{*}} of the correlation matrix according to 
  
  \deqn{R^{*} = \lambda T +  (1-\lambda) R }
  
  where \eqn{R} is the usual empirical correlation matrix and the target \eqn{T} is the unit diagonal
  matrix.  The shrinkage intensity \eqn{\lambda^{*}} for which the MSE of \eqn{R^{*}} is minimal
  is estimated by
  
   \deqn{\lambda^{*} = \sum_{i \neq j} Var(r_{ij}) / \sum_{i \neq j} r_{ij}^2  .}
 
  This is a special case of the analytic formula by Ledoit and Wolf (2003)
  for the optimal shrinkage (note that in contrast to Ledoit-Wolf here it is applied
  to the correlation rather than the covariance matrix).

  \code{cov.shrink} computes the corresponding full covariance matrix
  on the basis of the shrunken correlation matrix  and the empirical variances.
  
   \code{pcor.shrink} computes partial correlations from \code{cor.shrink}
  via \code{\link{cor2pcor}}. 
  
      
  These shrinkage estimator are especially useful in a 'small n, large p' setting - this is often
  encountered, e.g., in genomics. For an extensive discussion please see Schaefer and Strimmer (2005).
  
}
\value{
  \code{cov.shrink} returns a covariance matrix.
  
  \code{cor.shrink} returns the corresponding correlation matrix. 
  
   \code{pcor.shrink} returns the partical correlation matrix. 
}
\author{
  Juliane Schaefer (\url{http://www.statistik.lmu.de/~schaefer/}) 
  and Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
  
}
\references{

  Ledoit, O., and Wolf. M. (2003).  Improved estimation of the covariance
  matrix of stock returns with an application to portfolio selection.
  \emph{J. Emp. Finance} \bold{10}:503-621. 

  Schaefer, J., and Strimmer, K. (2005).  A shrinkage approach to large-scale
  covariance estimation and implications for functional genomics. 
  Statist. Appl. Genet. Mol. Biol.\bold{4}:32.
  (\url{http://www.bepress.com/sagmb/vol4/iss1/art32/})
  
  
  
}
\seealso{\code{\link{cor2pcor}}}

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
s1 <- cov(X)
s2 <- cov.shrink(X)

# squared error
sum((s1-sigma)^2)
sum((s2-sigma)^2)

# varcov produces the same results as cov
vc <- varcov(X)
sum(abs(vc$S-s1))

# compare positive definiteness
is.positive.definite(s1)
is.positive.definite(s2)
is.positive.definite(sigma)

# compare ranks and condition
rank.condition(s1)
rank.condition(s2)
rank.condition(sigma)

# compare eigenvalues
e1 <- eigen(s1, symmetric=TRUE)$values
e2 <- eigen(s2, symmetric=TRUE)$values
e3 <- eigen(sigma, symmetric=TRUE)$values
m <-max(e1, e2, e3)
yl <- c(0, m)

par(mfrow=c(1,3))
plot(e1,  main="empirical")
plot(e2,  ylim=yl, main="shrinkage")
plot(e3,  ylim=yl, main="true")
par(mfrow=c(1,1))

}
\keyword{multivariate}
