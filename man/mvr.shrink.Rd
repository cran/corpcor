\name{mvr.shrink}
\alias{mvr.shrink}
\alias{mvr.predict}


\title{Multivariate Shrinkage Regression}
\description{
  \code{mvr.shrink} implements an analytic shrinkage approach 
  fro multivariate linear regression that can be applied
  for "small n, large p" data.
  
  \code{mvr.predict} is a function to predict the response given 
  the (estimated) regression coefficients.
}

\usage{
mvr.shrink(x, y, lambda, lambda.var, w, verbose=TRUE)
mvr.predict(coef, x)
}

\arguments{
  \item{x}{a predictor matrix.}
  \item{y}{a response vector or matrix.  If \code{y} is missing then
           each individual column of \code{x} serves as
	    univariate response vector that is regressed in turn 
	    against the remaining predictors.}
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
            (\code{w = rep(1/n, n)} with \code{n = nrow(x)}).}	 
  \item{verbose}{output some status messages while computing (default: TRUE)}
  \item{coef}{the matrix of regression coefficients, as output by \code{mvr.shrink()}.}
  
}
\details{
   
   Most approaches to regularized regression focus solely on improving the estimated 
   covariance among the predictors (this is for example the case for ridge
   regression, lasso regression, elastic net, and its Bayesian cousins).
   
   In contrast, the shrinkage approach implemented in \code{mvr.shrink()} treats
   all response and predictor variables on an equal basis and regularizes the 
   \emph{joint} covariance matrix. Furthermore, if one considers both the 
   response and the predictors as random variables and estimates the best linear
   predictor with minimum MSE, then one arrives at a simple relationship
   between the shrinkage regression coefficients and shrunken partial covariance
   matrix.  For the specific details please see Opgen-Rhein and Strimmer (2006).
   
   In simulations this approach to multivariate shrinkage regression is generally
   more efficient than  standard least-squares and other regularized regression
   methods such as ridge regression.  
}
\value{
  \code{mvr.shrink} returns a matrix with the estimated regression coefficients.

  \code{mvr.predict} returns a response matrix.
}
\author{
  Rainer Opgen-Rhein (\url{http://www.statistik.lmu.de/~opgen/})
  and Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
  
}
\references{

  Opgen-Rhein, R., and Strimmer, K. (2006).  Beyond Gauss-Markov and ridge:
  an effective approach for high-dimensional regression and causal inference.
  Submitted.

  Schaefer, J., and Strimmer, K. (2005).  A shrinkage approach to large-scale
  covariance estimation and implications for functional genomics. 
  Statist. Appl. Genet. Mol. Biol.\bold{4}:32.
  (\url{http://www.bepress.com/sagmb/vol4/iss1/art32/})
   
}
\seealso{\code{\link{cov.shrink}}}

\examples{
# load corpcor library
library("corpcor")


# example data
data(longley) 
x     <- as.matrix(longley[ 1:13, -2])
y     <- as.matrix(longley[ 1:13,  2, drop=FALSE]) 
xtest <- as.matrix(longley[14:16, -2])
ytest <- as.matrix(longley[14:16,  2, drop=FALSE]) 


# least squares regression is recovered if
# shrinkage intensities are all set to zero
lm(y ~ x)
coefs.ols <- mvr.shrink(x, y, lambda=0, lambda.var=0)
coefs.ols


# shrinkage regression coefficients 
# note that these are quite different!
coefs.shrink <- mvr.shrink(x, y)
coefs.shrink


# prediction
mvr.predict(coefs.ols, xtest)
mvr.predict(coefs.shrink, xtest)


# if no response matrix is given, then each predictor 
# is regressed in turn against all others
coefs <- mvr.shrink(longley)
coefs


# connection between partial correlations and regression coefficients
# (note that the mvr.shrink() function is specifically constructed so that 
# this relationship holds for the shrinkage estimates)

beta <- coefs[,-1]    # shrinkage regression coefficients w/o intercept
pc <- pcor.shrink(longley) # shrunken partial correlations
k <- 3
l <- 4
# partial correlations as computed by pcor.shrink()
pc[k,l] 
# partial correlations obtained from mvr.shrink() regression coefficients 
sqrt(beta[k,l] * beta[l,k]) * sign(beta[l,k]) 

}
\keyword{multivariate}
