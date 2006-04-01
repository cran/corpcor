\name{smtools}
\alias{sm2vec}
\alias{sm.indexes}
\alias{vec2sm}

\title{Some Tools for Symmetric Matrices}

\usage{
sm2vec(m, diag = FALSE)
sm.indexes(m, diag = FALSE)
vec2sm(vec, diag = FALSE, order = NULL)
}

\arguments{
  \item{m}{symmetric matrix}
  \item{diag}{logical. Should the diagonal be included?}
  \item{vec}{vector of unique elements from a symmetric matrix}
  \item{order}{order of the entries in \code{vec}}
}
\description{
  \code{sm2vec} takes a symmetric matrix and puts
  the lower triagonal entries into a vector (cf. \code{\link{lower.tri}}).

  \code{sm.indexes} gives the corresponding x-y-indexes for each entry
  in the vector produced by \code{sm2vec}.
  
  \code{vec2sm} reverses the operation by  \code{sm2vec} and turns the
  vector back in a symmetric matrix. If \code{diag=FALSE} the
  diagonal of the resulting matrix will consist of NAs.  If \code{order}
  is supplied then the input vector \code{vec} will first be rearranged accordingly.
 }

\value{
  A vector (\code{sm2vec}), a two-column matrix with indexes (\code{sm.indexes}),
  or a symmetric matrix (\code{vec2sm}).
}


\author{
  Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
}

\seealso{\code{\link{lower.tri}}.}

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

# convert into vector (including diagonals
v <- sm2vec(m.cov, diag=TRUE)
v.idx <- sm.indexes(m.cov, diag=TRUE)
v
v.idx

# put back to symmetric matrix
vec2sm(v, diag=TRUE)

# vector not in the original order
sv <- sort(v)
sv
ov <- order(v)
ov
vec2sm(sv, diag=TRUE, order=ov)
}
\keyword{utilities}
