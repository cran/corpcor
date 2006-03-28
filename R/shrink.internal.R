### shrink.internal.R  (2006-03-27)
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

# function to compute shrinkage variance vector
#  - wm: weighted moments of original x matrix
#  - xc: *centered* data matrix, 
#  - lambda.var = 0: don't shrink
#    lambda.var > 0: shrink with given lambda
#    lambda.var < 0: shrink with estimated lambda  
#  - w:  data weights

pvt.svar <- function(wm, xc, lambda.var, w, verbose=TRUE)
{  
  z <- pvt.get.lambda(xc, lambda.var, w, verbose=verbose, type="variance")
        
  # shrunken variances
  sv <- z$lambda.var*mean(wm$var) + (1-z$lambda.var)*wm$var
       
  attr(sv, "lambda.var") <- z$lambda.var
  attr(sv, "lambda.var.estimated") <- z$lambda.var.estimated
  
  return(sv)   
}    



# function to compute shrinkage correlation matrix 
#  - xs: scaled data matrix, 
#  - lambda = 0: don't shrink
#    lambda > 0: shrink with given lambda
#    lambda < 0: shrink with estimated lambda  
#  - w:  data weights

pvt.scor <- function(xs, lambda, w, verbose=TRUE)
{  
  z <- pvt.get.lambda(xs, lambda, w, verbose=verbose, type="correlation")
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
    }
    
    # set all diagonal entries to 1
    diag(r) <- 1 
    
  }  
  attr(r, "lambda") <- z$lambda
  attr(r, "lambda.estimated") <- z$lambda.estimated
 
  return(r)   
}    

# compute the *inverse* of the shrinkage correlation matrix

# this procedure exploits the Woodbury identity to efficiently
# compute the inverse of the correlation shrinkage estimator
# directly 

pvt.invscor <- function(wm, xs, lambda, w, verbose=TRUE, type="correlation")
{
  z <- pvt.get.lambda(xs, lambda, w, verbose=verbose)
  p <- dim(xs)[2]
    
  if (z$lambda == 1)
  {
    invr <- diag(p)
  }
  else
  {
    # number of zero-variance variables
    zeros <- (wm$var==0.0)
    
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
      if (m < p-sum(zeros)) 
        warning(paste("Estimated correlation matrix doesn't have full rank",
          "- pseudoinverse used for inversion."), call. = FALSE)    
       
      invr <- svdxs$v %*% iC %*% t(svdxs$v)
    }
    else # use a special case of Woodbury matrix identity for inversion
    {
      F <- solve(z$lambda*iC + diag(m))
      invr <- (diag(p) - svdxs$v %*% F %*% t(svdxs$v))/z$lambda 
    }
 
    # set all diagonal entries corresponding to zero-variance variables to 1
    diag(invr)[zeros] <- 1
 
  }  
  
  attr(invr, "lambda") <- z$lambda
  attr(invr, "lambda.estimated") <- z$lambda.estimated
  
  return( invr )
}



# returns lambda to be used for shrinkage
pvt.get.lambda <- function(x, lambda, w, verbose=TRUE, type=c("correlation", "variance"))
{
  # if lambda = 0: don't shrink
  # if lambda > 0: shrink with given lambda
  # if lambda < 0: shrink with estimated lambda  (the default)
 
  type <- match.arg(type)
  
  if (type == "correlation")
  {
     kind <- "correlation matrix"
     func <- "C_corlambda"
     
     # note: x needs to be the *scaled* data matrix
  }
  
  if (type == "variance")
  {
     kind <- "variance vector"
     func <- "C_varlambda"
     
     # note: x needs to be the *centered* data matrix
  }
  
    
  if (lambda < 0)
  { 
    if (verbose)
    {
      cat(paste("Determining optimal shrinkage intensity (", kind, ") ...\n", sep=""))     
    }
    
    # estimate optimal shrinkage intensity 
    # target: correlations/covariances -> 0  
    lambda <- .C(func,
            as.double(x),
	    as.integer( dim(x)[1] ),
	    as.integer( dim(x)[2] ),
	    as.double(w),
	    lambda=double(1), PACKAGE="corpcor", DUP=FALSE)$lambda

    lambda.estimated <- TRUE
      
    if (verbose)
    {
      cat(paste("Estimated shrinkage intensity (", kind, "): ", round(lambda, 4), "\n", sep=""))     
    }
  }
  else
  {
    if (lambda > 1) lambda <- 1
    lambda.estimated <- FALSE
    
    if (verbose)
    {
      cat(paste("Specified shrinkage intensity (", kind, "): ", round(lambda, 4), "\n", sep=""))     
    }    
  }
  
  if (type == "correlation") 
  {
     return (list(lambda=lambda, lambda.estimated=lambda.estimated))
  }


  if (type == "variance") 
  {
     return (list(lambda.var=lambda, lambda.var.estimated=lambda.estimated))
  }
	       
}




