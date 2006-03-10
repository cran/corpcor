### shrink.internal.R  (2006-03-09)
###
###    Non-public functions used in the covariance shrinkage estimator 
###    
###
### Copyright 2005-06 Korbinian Strimmer
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


# make sure we have weights that sum up to one
pvt.check.w <- function(w, n)
{
   if (missing(w))
   {
     w <- rep(1/n, n)   # equal weights
   }
   else
   {
     if (length(w) != n)
     {
       stop("Weight vector has incompatible length", call. = FALSE)
     } 
     
     sw <- sum(w)
     if (sw != 1) 
       w <- w/sw      # make weights sum up to 1
   } 
 
   return(w)
}



# function to compute shrinkage correlation matrix 
#  - xs: scaled data matrix, 
#  - lambda = 0: don't shrink
#    lambda > 0: shrink with given lambda
#    lambda < 0: shrink with estimated lambda  
#  - w:  data weights

pvt.scor <- function(xs, lambda, w, verbose=TRUE)
{  
  z <- pvt.get.lambda(xs, lambda, w, verbose=verbose)
  if (z$lambda == 1)
  {
    p <- dim(xs)[2]
    r <- diag(p)
  }
  else
  {
    # bias correction factor
    h1 <- 1/(1-sum(w*w))   # for w=1/n this equals the usual h1=n/(n-1)
 
    # unbiased empirical estimator
    # for w=1/n  the following  would simplify to:  r <- 1/(n-1)*crossprod(xs)
    #r <- h1 * t(xs) %*% diag(w) %*% xs
    r <- h1 * t(xs) %*% sweep(xs, 1, w, "*") # sweep requires less memory
   
    if (z$lambda > 0) 
    {
      # shrink off-diagonal elements
      r <- (1-z$lambda)*r
      diag(r) <- 1
    }
  }  
  attr(r, "overshrinkage") <- z$overshrinkage
  attr(r, "lambda") <- z$lambda
  attr(r, "lambda.estimated") <- z$lambda.estimated
 
  return(r)   
}    

# compute the *inverse* of the shrinkage correlation matrix

# this procedure exploits the Woodbury identity to efficiently
# compute the inverse of the correlation shrinkage estimator
# directly 

pvt.invscor <- function(xs, lambda, w, verbose=TRUE)
{
  z <- pvt.get.lambda(xs, lambda, w, verbose=verbose)
  p <- dim(xs)[2]
    
  if (z$lambda == 1)
  {
    invr <- diag(p)
  }
  else
  {
    svdxs <- fast.svd(xs)
    m <- length(svdxs$d)  # rank of xs
           
    # bias correction factor
    h1 <- 1/(1-sum(w*w))   # for w=1/n this equals the usual h1=n/(n-1)
    
    UTWU <- t(svdxs$u) %*% sweep(svdxs$u, 1, w, "*") #  t(U) %*% diag(w) %*% U
    C <- sweep(sweep(UTWU, 1, svdxs$d, "*"), 2, svdxs$d, "*") # D %*% UTWU %*% D
    C <- (1-z$lambda) * h1 * C
    iC <- solve(C)
    
    # note: C is of size m x m, and diagonal if w=1/n
     
    if (lambda==0.0) # use eigenvalue decomposition for inversion
    {
      if (m < p) warning(paste("Estimated correlation matrix doesn't have full rank",
      "- pseudoinverse used for inversion."), call. = FALSE)    
       
      invr <- svdxs$v %*% iC %*% t(svdxs$v)
    }
    else # use a special case of Woodbury matrix identity for inversion
    {
      F <- solve(z$lambda*iC + diag(m))
      invr <- (diag(p) - svdxs$v %*% F %*% t(svdxs$v))/z$lambda 
    }
  }  
  
  attr(invr, "overshrinkage") <- z$overshrinkage
  attr(invr, "lambda") <- z$lambda
  attr(invr, "lambda.estimated") <- z$lambda.estimated
  
  return( invr )
}



# returns lambda to be used for shrinkage
pvt.get.lambda <- function(xs, lambda, w, verbose=TRUE)
{
  # if lambda = 0: don't shrink
  # if lambda > 0: shrink with given lambda
  # if lambda < 0: shrink with estimated lambda  (the default)
 
  if (lambda < 0)
  { 
    lambda <- pvt.estimate.lambda(xs, w)
    lambda.estimated <- TRUE
      
    if (verbose)
    {
      cat(paste("Estimated shrinkage intensity lambda: ", round(lambda, 4), "\n"))     
    }
  }
  else
  {
    lambda.estimated <- FALSE
  }
  
  if (lambda > 1)
  {  
    overshrinkage <- TRUE
    warning(paste("Overshrinking: intensity lambda set to 1", 
     "(allowed range: 0-1)"), call. = FALSE, immediate. = FALSE)
    lambda <- 1
  }
  else
    overshrinkage <- FALSE
      
  return (list(lambda=lambda, 
               lambda.estimated=lambda.estimated,
	       overshrinkage=overshrinkage
	       ))
}




######################################################################
# this internal function may profit from being rewritten in C ...
######################################################################

# input:  scaled data matrix  (standardization is NOT checked)
#         weights of each data point
#         estimated lambda 
#
# Note: this procedure does NOT store nor generate 
#       the complete pxp correlation matrix
pvt.estimate.lambda <- function(xs, w)
{
  # bias correction factors
  h1 <- 1/(1-sum(w*w))   # for w=1/n this equals the usual h1=n/(n-1)
  h2 <- h1*h1*sum(w*w)   # for w=1/n this equals h2=n/(n-1)^2 
  h3 <- h1*h2            # for w=1/n this equals h3=n^2/(n-1)^3 
 
  p <- dim(xs)[2]
   
  offdiagsum.rij.2 <- 0
  offdiagsum.v.rij <- 0

  for (i in 1:(p-1))
  {
    for (j in (i+1):p)
    {
       xsij <- xs[,i]*xs[,j]       
       rr <- h1*sum(w*xsij)
       xsijc <- xsij-rr
       vr <- h3*sum(w*xsijc*xsijc)
                      
       #rr <- h1*weighted.mean(xsij, w)
       #vr <- h2*weighted.var(xsij, w)
		      		      
       offdiagsum.v.rij <- offdiagsum.v.rij + vr
       offdiagsum.rij.2 <- offdiagsum.rij.2 + rr*rr
    }
  }

  lambda <- offdiagsum.v.rij/offdiagsum.rij.2

  return (lambda)
}




