\name{fast.svd}
\alias{fast.svd}

\title{Fast Singular Value Decomposition}
\usage{
fast.svd(m, tol)
}
\arguments{
  \item{m}{matrix}
  \item{tol}{tolerance - singular values larger than
             tol are considered non-zero (default value:
	     \code{tol = max(dim(m))*max(D)*.Machine$double.eps})
	     }
}
\description{
  Any rectangular real matrix M can be decomposed as

                             \deqn{M = U D V^{'},}{M = U D V',}

  where U and V are orthogonal, V' means V transposed, and 
  D is a diagonal matrix with the singular values (see \code{\link{svd}}).
  
  If the matrix M is either "fat" (small n, large p) or "thin" (large n,
  small p) then then decomposition of M into the U, D, V matrices can be
  greatly sped up by computing the SVD of either M M' (fat matrices) or
  M' M (thin matrices), rather than that of M.
  
  This is the trick that allows \code{fast.svd} to be substantially faster 
  than the native \code{\link{svd}} for fat and thin matrices.  For
  near-square matrices \code{fast.svd} simply uses the standard \code{\link{svd}}
  function.
  
  In contrast to the native \code{\link{svd}} function note that
  \code{fast.svd} returns only *positive* singular values
  (i.e. rank D equals rank M) and their associated singular vectors. 
   Also note that the latter
  may differ in sign from  \code{\link{svd}}. 
}

\details{
  For "fat" M (small n, large p) the SVD decomposition of M M' yields
  
            \deqn{M M^{'} = U D^2 U}{M M' = U D^2 U'}
	     
  As the matrix M M' has dimension n x n only, this is faster to compute
  than SVD of M.  The V matrix is subsequently obtained by
  
             \deqn{V = M^{'} U D^{-1}}{V = M' U D^(-1)}  
	            
  Similarly, for "thin" M (large n, small p), the decomposition of M' M
  yields
  
             \deqn{M^{'} M = V D^2 V^{'}}{M' M = V D^2 V'}
	    
  which is also quick to compute as M' M has only dimension p x p.  The 
  U matrix is then computed via
  
             \deqn{U = M V D^{-1}}{U = M V D^(-1)}
}

\value{
  A list with the following components:
  \item{d}{a vector containing the \emph{positive} singular values }
  \item{u}{a matrix with the corresponding left singular vectors}
  \item{v}{a matrix with the corresponding left singular vectors}
}


\author{
  Korbinian Strimmer (\url{http://www.statistik.lmu.de/~strimmer/}).
}
\seealso{\code{\link{svd}}, \code{\link{solve}}.}

\examples{
# load corpcor library
library("corpcor")


# generate a "fat" data matrix
n <- 50
p <- 5000
X <- matrix(rnorm(n*p), n, p)

# compute SVD
system.time( s1 <- svd(X) ) 
system.time( s2 <- fast.svd(X) )


eps <- 1e-10
sum(abs(s1$d-s2$d) > eps)
sum(abs(abs(s1$u)-abs(s2$u)) > eps)
sum(abs(abs(s1$v)-abs(s2$v)) > eps)
}
\keyword{algebra}
