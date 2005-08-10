### fast.svd.R  (2005-07-19)
###
###    Efficient Computation of the Singular Value Decomposition
###
### Copyright 2003-05 Korbinian Strimmer
###
###
### This file is part of the `corpcor' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

# private functions

# this works just like svd, with the difference that only LAPACK
# methods are used, and that there is an automatic fallback
# to algorithm "dgesvd" if algorithm "dgesdd" throws an error.

# note that LAPACK.svd returns "v" not "vt" as La.svd


LAPACK.svd <- function(x, nu = min(n, p), nv = min(n, p))
{
    if (!is.numeric(x) && !is.complex(x))
        stop("argument to 'LAPACK.svd' must be numeric or complex")
    if (any(!is.finite(x)))
        stop("infinite or missing values in 'x'")
         
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if (!n || !p)
        stop("0 extent dimensions")


    if( is.complex(x) ) # if complex we have to use the "dgesvd" algorithm
    {
        method = "dgesvd"
        res <- La.svd(x, nu=nu, nv=nv, method=method) 
	
    }
    else # but otherwise we try the "dgesdd" algorithm first (10x faster)
    {
        method = "dgesdd"
        res <- try( La.svd(x, nu=nu, nv=nv, method=method), silent=TRUE )
  
        if ( class(res) == "try-error" || any(res$d < 0) ) 
	# and use the "dgesvd" algorithm if there is an error,
	# or if there are negative singular values (happens on AMD64/ATLAS!!)  
        {
            method = "dgesvd"
            res <- La.svd(x, nu=nu, nv=nv, method=method)
        }
    }

    if  (is.complex(x))
    {
         out <- list(d = res$d, u = if (nu) res$u, v = if (nv) Conj(t(res$vt)))
    }
    else
    {
         out <- list(d = res$d, u = if (nu) res$u, v = if (nv) t(res$vt))
    }

    attr(out, "svd.algorithm") <- method
    return(out)
}



# standard svd returning only positive singular values
positive.svd <- function(m, tol)
{
  s <- LAPACK.svd(m)
  
  if( missing(tol) ) 
      tol <- max(dim(m))*max(s$d)*.Machine$double.eps
  Positive <- s$d > tol

  return(list(
      d=s$d[Positive],
      u=s$u[, Positive, drop=FALSE],
      v=s$v[, Positive, drop=FALSE]
      ))
}

# fast computation of svd(m) if n << p  
# (n are the rows, p are columns)
nsmall.svd <- function(m, tol)
{
   B <- m %*% t(m)     # nxn matrix
   s <- LAPACK.svd(B,nv=0)    # of which svd is easy..

   # determine rank of B  (= rank of m)
   if( missing(tol) ) 
      tol <- dim(B)[1]*max(s$d)*.Machine$double.eps 
   Positive <- s$d > tol                            
           
   # positive singular values of m  
   d <- sqrt(s$d[Positive])
      
   # corresponding orthogonal basis vectors
   u <- s$u[, Positive, drop=FALSE]
   v <- crossprod(m, u) %*% diag(1/d, nrow=length(d))   
  
   return(list(d=d,u=u,v=v))
}

# fast computation of svd(m) if n >> p  
# (n are the rows, p are columns)
psmall.svd <- function(m, tol)
{
   B <- crossprod(m)   # pxp matrix
   s <- LAPACK.svd(B,nu=0)    # of which svd is easy..

   # determine rank of B  (= rank of m)
   if( missing(tol) ) 
      tol <- dim(B)[1]*max(s$d)*.Machine$double.eps 
   Positive <- s$d > tol                            
           
   # positive singular values of m  
   d <- sqrt(s$d[Positive])
      
   # corresponding orthogonal basis vectors
   v <- s$v[, Positive, drop=FALSE]
   u <- m %*% v %*% diag(1/d, nrow=length(d))
  
   return(list(d=d,u=u,v=v))
}


# public functions

# fast computation of svd(m)

# note that the signs of the columns vectors in u and v
# may be different from that given by svd()

# note that also only positive singular values are returned

fast.svd <- function(m, tol)
{  
  n <- dim(m)[1]
  p <- dim(m)[2]
 
 
  EDGE.RATIO <- 2 # use standard SVD if matrix almost square
  if (n > EDGE.RATIO*p)
  {
     return(psmall.svd(m,tol))
  }
  else if (EDGE.RATIO*n < p)
  {  
     return(nsmall.svd(m,tol)) 
  }
  else # if p and n are approximately the same
  {
     return(positive.svd(m, tol))
  }
}


