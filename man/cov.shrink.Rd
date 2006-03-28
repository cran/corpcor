\name{cov.shrink}
\alias{cov.shrink}
\alias{cor.shrink}
\alias{var.shrink}


\title{Shrinkage Estimates of Covariance and Correlation}
\description{
  The functions \code{cov.shrink} and \code{cor.shrink} compute shrinkage estimates of
  covariance and correlation, as as suggested in Schaefer and Strimmer (2005).
    
  In comparison with the standard empirical estimates 
  (\code{\link{cov}} and \code{\link{cor}}) the shrinkage estimates have a number
  of advantages
  \enumerate{
    \item are always  positive definite,
    \item  well conditioned (so that the inverse always exists), and 
    \item exhibit (sometimes dramatically) better mean squared error.}
  Furthermore, they are inexpensive to compute and do not require any
  tuning parameters (the shrinkage intensity is analytically estimated from the data).
  
  As an extra benefit, the shrinkage estimators have a form that can be *very* efficiently inverted
  using the Woodbury matrix identity, even if the number of variables is large.  
  This approach is implemented in the functions 
  \\code{\link{invcov.shrink}} and \code{\link{invcor.shrink}} and is much faster than directly inverting
  the matrix output by \code{cov.shrink} and \code{cor.shrink}, respectively.
}
\usage{
var.shrink(x, lambda.var, w, verbose=TRUE)
cor.shrink(x, lambda, w, verbose=TRUE)
cov.shrink(x, lambda, lambda.var, w, verbose=TRUE)
}
\arguments{
  \item{x}{a data matrix}
  \item{lambda}{the correlation shrinkage intensity (range 0-1).
                If \code{lambda} is not specified (the default) it is estimated
		using an analytic formula from Schaefer and Strimmer (2005) 
		- see details below.  
		For \code{lambda=0} the empirical correlations are recovered.}
  \item{lambda.var}{the variance shrinkage intensity (range 0-1). 
                If \code{lambda.var} is not specified (the default) it is estimated
	        using an analytic formula from Schaefer and Strimmer (2005) 
                - see details below.  
	        For \code{lambda.var=0} the empirical variances are recovered.}	 
  \item{w}{optional: weights for each data point - if not specified uniform weights are assumed
            (\code{w = rep(1/n, n)} with \code{n = dim(x)[1]}).}	 
  \item{verbose}{output some status messages while computing (default: TRUE)}
}
\details{
   
   \code{var.shrink} computes the empirical variances of each variables, and shrinks them towards
   the average of all variances.  The shrinkage intensity is estimated using Eq. 6 from 
   Schaefer and Strimmer (2005), which leads to
   
   
      \deqn{\lambda_{var}^{*} = ( (1-1/p) \sum_{k=1}^p Var(s_{kk}) )/ \sum_{k=1}^p (s_{kk} - \bar{s})^2  .}

  Similarly   \code{cor.shrink} computes a shrinkage estimate of the correlation matrix by
  shrinking the empirical correlations towards  the identity matrix.  
  In this case the shrinkage intensity is
  
  
      \deqn{\lambda^{*} = \sum_{k \neq k} Var(r_{kl}) / \sum_{k \neq l} r_{kl}^2  .}
 

  \code{cov.shrink} computes the corresponding full covariance matrix
  on the basis of the shrunken correlation matrix  and the shrunken variances.
    
  These shrinkage estimators are especially useful in a 'small n, large p' setting - this is often
  encountered, e.g., in genomics. For an extensive discussion please see Schaefer and Strimmer (2005).
}
\value{
  \code{cov.shrink} returns a covariance matrix.
  
  \code{cor.shrink} returns the corresponding correlation matrix. 
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
\seealso{\code{\link{invcov.shrink}}, \code{\link{pcor.shrink}}, \code{\link{cor2pcor}}}

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
