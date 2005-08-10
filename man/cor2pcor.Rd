\name{cor2pcor}
\alias{cor2pcor}
\alias{pcor2cor}

\title{Compute Partial Correlation from Correlation Matrix -- and Vice Versa}
\usage{
cor2pcor(m, exact.inversion=TRUE, check.eigenvalues=TRUE, tol)
pcor2cor(m, exact.inversion=TRUE, check.eigenvalues=TRUE, tol)
}
\arguments{
  \item{m}{covariance matrix or (partial) correlation matrix}
  \item{check.eigenvalues}{if TRUE the input matrix is checked for positive definiteness
        (i.e. whether all eigenvalues are strictly postive)}
  \item{exact.inversion}{determines whether the inverse is computed
        exactly (using \code{\link{solve}}) or via \code{\link{pseudoinverse}}  }
   
  \item{tol}{tolerance - singular values larger than
             tol are considered non-zero (default value:
	     \code{tol = max(dim(m))*max(D)*.Machine$double.eps}).
	     This values is used both for  testing positive definiteness
	     as well as for computing the pseudoonverse.
	     }
}
\description{
  \code{cor2pcor} computes the pairwise 
  \emph{partial} correlation coefficients from either a correlation 
  or a covariance matrix. 

  \code{pcor2cor} takes a partial correlation matrix and computes
  the corresponding correlation matrix.
  
}

\details{
  The partial correlations represent the direct
  interactions between two variables, with the indirect effects of all
  remaining variables removed.
  
  For computing partial correlation coefficients from data use the
  function \code{\link{pcor.shrink}}.
}

\value{
  A matrix with the pairwise partial correlation coefficients
  (\code{cor2pcor} and \code{pcor}) or with pairwise
  correlations (\code{pcor2cor})
}


\author{
  Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
}
\references{
  Whittaker J. (1990).  Graphical Models in Applied Multivariate Statistics.
   John Wiley, Chichester.
}

\seealso{\code{\link{cor}}, \code{\link{pseudoinverse}}.}

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

# corresponding correlation matrix
m.cor.1 <- cov2cor(m.cov)
m.cor.1

# compute partial correlations (from covariance matrix)
m.pcor.1 <- cor2pcor(m.cov)
m.pcor.1

# compute partial correlations (from correlation matrix)
m.pcor.2 <- cor2pcor(m.cor.1)
m.pcor.2

zapsmall( m.pcor.1 ) == zapsmall( m.pcor.2 )

# backtransformation
m.cor.2 <- pcor2cor(m.pcor.1)
m.cor.2
zapsmall( m.cor.1 ) == zapsmall( m.cor.2 )

}
\keyword{multivariate}
