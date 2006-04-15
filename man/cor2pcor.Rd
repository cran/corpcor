\name{cor2pcor}
\alias{cor2pcor}
\alias{pcor2cor}
\alias{cov2pcov}
\alias{pcov2cov}

\title{Compute Partial Correlation from Correlation Matrix -- and Vice Versa}
\usage{
cor2pcor(m, tol)
pcor2cor(m, tol)
cov2pcov(m, tol)
pcov2cov(m, tol)
}
\arguments{
  \item{m}{covariance matrix or (partial) correlation matrix}   
  \item{tol}{tolerance - singular values larger than
             tol are considered non-zero (default value:
	     \code{tol = max(dim(m))*max(D)*.Machine$double.eps}).
	     This parameter is needed for  the singular
	     value decomposition on which \code{\link{pseudoinverse}} is based.
	     }
}
\description{
  \code{cor2pcor} computes the pairwise 
  \emph{partial} correlation coefficients from either a correlation 
  or a covariance matrix. 

  \code{pcor2cor} takes either a partial correlation matrix or 
   a partial covariance matrix as input,
   and computes from it the corresponding correlation matrix.
   
   \code{cov2pcov} transform a covariance matrix into a partial
   covariance matrix, and  \code{pcov2cov} is the corresponding
   reciprocal operation.
}

\details{
  The partial covariance is  the inverse of the covariance matrix,
  with the signs of the off-diagonal elements reversed.  The partial
  correlation is the standardized partial covariance matrix.

  In graphical Gaussian models the partial correlations represent the 
  direct interactions between two variables, conditioned on all
  remaining variables.
  
  In all the above functions the \code{\link{pseudoinverse}} is employed
  for inversion - hence even singular covariances (with some
  zero eigenvalues) may be used.  However, a better option may be to
  estimate a positive definite covariance matrix using 
  \code{\link{cov.shrink}}.
  
  Note that for efficient computation of partial correlation coefficients from some 
  data  \code{x} use  \code{pcor.shrink(x)}  and \emph{not} 
  \code{cor2pcor(cor.shrink(x))}.
}

\value{
  A matrix with the pairwise partial correlation coefficients
  (\code{cor2pcor}), with pairwise
  correlations (\code{pcor2cor}), with the partial covariances (\code{cov2pcor}),
  or the covariances (\code{pcov2cov}).
}


\author{
  Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
}
\references{
  Whittaker J. (1990).  Graphical Models in Applied Multivariate Statistics.
   John Wiley, Chichester.
}

\seealso{\code{\link{pcor.shrink}}, \code{\link{pseudoinverse}}.}

\examples{
# load corpcor library
library("corpcor")

# covariance matrix
m.cov <- rbind(
 c(3,1,1,0),
 c(1,3,0,1),
 c(1,0,2,0),
 c(0,1,0,2)
)
m.cov


# partial covariance matrix
m.pcov <- cov2pcov(m.cov)
m.pcov


# corresponding correlation matrix
m.cor.1 <- cov2cor(m.cov)
m.cor.1

# compute partial correlations (from covariance matrix)
m.pcor.1 <- cor2pcor(m.cov)
m.pcor.1

# compute partial correlations (from correlation matrix)
m.pcor.2 <- cor2pcor(m.cor.1)
m.pcor.2

# compute partial correlations (from partial covariance matrix)
m.pcor.3 <- cov2cor(m.pcov)
m.pcor.3


zapsmall( m.pcor.1 ) == zapsmall( m.pcor.2 )

# backtransformation
m.cor.2 <- pcor2cor(m.pcor.1)
m.cor.2
zapsmall( m.cor.1 ) == zapsmall( m.cor.2 )

}
\keyword{multivariate}
