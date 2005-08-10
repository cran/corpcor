\name{cov.shrink}
\alias{cov.shrink}
\alias{cor.shrink}
\alias{pcor.shrink}

\title{Shrinkage Estimates of Covariance and (Partial) Correlation}
\description{
  The functions \code{cov.shrink}, \code{cor.shrink}, and \code{pcor.shrink} implement a shrinkage approach 
  to estimate covariance and (partial) correlation matrices.
  
  The advantages of using this approach in comparison with the standard empirical estimates 
  (\code{\link{cov}} and \code{\link{cor}}) are that 
   1) the shrinkage estimates are always  positive definite, 2)  well conditioned (so that the inverse
  always exists), and 3) exhibit (sometimes dramatically) better mean squared error.
  Furthermore, they are inexpensive to compute and do not need the specification of any
  parameters (the shrinkage intensity is analytically estimated from the data).
  
}
\usage{
cov.shrink(x, lambda, verbose=TRUE)
cor.shrink(x, lambda, verbose=TRUE)
pcor.shrink(x, lambda, verbose=TRUE)
}
\arguments{
  \item{x}{a data matrix}
  \item{lambda}{the shrinkage intensity (range 0-1). It is is not specified (the default) lambda is
         chosen such that the resulting shrinkage estimate has minimal MSE (see below for details).}
  \item{verbose}{report progress while computing (default: TRUE)}
}
\details{
   
  Let \eqn{S} be the usual empirical unbiased estimator of the true covariance matrix 
  \eqn{\Sigma} with \eqn{E(S)=\Sigma},
  and the target \eqn{T} an unbiased estimator of a structured covariance matrix 
  (with \eqn{E(T)=\Theta}.  Then a shrinkage
  estimator \eqn{S^{*}} can be constructed by setting  
  
      
     \deqn{S^{*} = \lambda T +  (1-\lambda) S .}
      
  
  If the shrinkage intensity  \eqn{\lambda} is large, then the resulting estimate 
  will be close the structured
  estimate \eqn{T}, otherwise the unstructured estimate \eqn{S} dominates.   
  
  
  Ledoit and Wolf (2003) suggest an analytic formula to estimate \eqn{\lambda} 
  from the observed data such that the mean squared error (MSE)
  of \eqn{S^{*}} is minimal. 
   
  
  In \code{cov.shrink} the target \eqn{T} is the diagonal matrix \eqn{diag(S)}. 
  The above equation then simplifies to a shrinkage estimator of the offdiagonal elements only, i.e.
  with \eqn{s^{*}_{ii} = s_{ii}}
  and \eqn{s^{*}_{ij} = (1-\lambda) s_{ij}} for \eqn{i \neq j}.
  It then can be shown (see Schaefer and Strimmer 2005) that an estimator for the optimal intensity is
  given by
  
  \deqn{\lambda^{*} = \sum_{i \neq j} Var(s_{ij}) / \sum_{i \neq j} s_{ij}^2  .}
    
    
  This shrinkage estimator are especially useful when computing sparse large-scale covariance
  matrices in a 'small n, large p' setting. These situations are often encountered
  in problems bioinformatics and statistical genomics (see Schaefer and Strimmer 2005 for examples).
}
\value{
  \code{cov.shrink} returns a matrix with the shrinkage estimate of the covariance matrix.
  
  \code{cor.shrink} returns the corresponding standardized matrix. 
  
   \code{pcor.shrink} returns the partical correlations computed from the output of \code{cov.shrink}. 
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
  covariance estimation and implications for functional genomics. Submitted to SAGMB. 
  
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
