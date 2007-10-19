### shrink.internal.R  (2007-10-19)
###
###    Non-public functions used in the covariance shrinkage estimator 
###    
###
### Copyright 2005-07 Korbinian Strimmer
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


# print function
print.shrinkage <- function(x, ...)
{
  attr(x, "class") <- NULL
  
  lambda <- attr(x, "lambda")
  lambda.estimated <- attr(x, "lambda.estimated")
  attr(x, "lambda") <- NULL 
  attr(x, "lambda.estimated") <- NULL 
  
  lambda.var <- attr(x, "lambda.var")
  lambda.var.estimated <- attr(x, "lambda.var.estimated")
  attr(x, "lambda.var") <- NULL 
  attr(x, "lambda.var.estimated") <- NULL 
  
  protect <- attr(x, "protect")
  attr(x, "protect") <- NULL

  spv <- attr(x, "spv")
  attr(x, "spv") <- NULL

  NextMethod("print", x, quote = FALSE, right = TRUE)
    
  cat("\n")
  if (!is.null(lambda.estimated))  
  {
    if (lambda.estimated) 
      le <- "(estimated)"
    else
      le <- "(specified)"  
    cat(paste("Shrinkage intensity lambda (correlation matrix):", round(lambda,4), le, "\n"))
  }
  
  if (!is.null(lambda.var.estimated))  
  {
    if (lambda.var.estimated)
      lve <- "(estimated)"
    else
      lve <- "(specified)"  
    cat(paste("Shrinkage intensity lambda.var (variance vector):", round(lambda.var, 4), lve, "\n"))
  }
  
  if (!is.null(protect))  
  {
     cat("Fraction of components with limited translation:", protect, "\n")
  }
 
  if (!is.null(spv))  
  {
     cat("Standardized partial variances (i.e. PVAR/VAR) are attached (attribute \"spv\").\n")
  }
}




# function to compute shrinkage variance vector
#  - x: data matrix, 
#  - lambda.var = 0: don't shrink
#    lambda.var > 0: shrink with given lambda
#    lambda.var < 0: shrink with estimated lambda  
#  - w:  data weights

pvt.svar <- function(x, lambda.var, w, protect=protect, verbose)
{    
  # center input matrix
  xs <- wt.scale(x, w, center=TRUE, scale=FALSE) 
  
  # compute variances 
  v <- wt.moments(xs, w)$var
  
  # compute target
  tgt <- median(v)
          
  # shrinkage estimate 
  z <- pvt.get.lambda(xs, lambda.var, w, 0, verbose=verbose, type="variance", tgt)   
  vs <- z$lambda.var*tgt + (1-z$lambda.var)*v

  if (protect > 0)
  {
    if (protect >= 1)
    {
      vs <- v
      protect=1
    }
    else
    {
      diff <- v-vs 
      adiff <- abs(diff)
      sdiff <- sign(diff)
      diff <- NULL
  
      M <- quantile(adiff, probs=c(1-protect))
  
      d <- adiff-M # soft thresholding
      d[d < 0] <- 0
      d <- d*sdiff
    
      vs <- vs + d  # add correction term
    }
  
    attr(vs, "protect") <- protect
  }

         
  attr(vs, "lambda.var") <- z$lambda.var
  attr(vs, "lambda.var.estimated") <- z$lambda.var.estimated
  attr(vs, "class") <- "shrinkage"
  
  return(vs)   
}    



# function to compute shrinkage correlation matrix 
#  - x: data matrix, 
#  - lambda = 0: don't shrink
#    lambda > 0: shrink with given lambda
#    lambda < 0: shrink with estimated lambda  
#  - w:  data weights

pvt.scor <- function(x, lambda, w, protect, verbose)
{  
  # standardize input matrix by standard deviations
  xs <- wt.scale(x, w, center=TRUE, scale=TRUE, scale.by="sd") 
  
  # bias correction factor
  h1 <- 1/(1-sum(w*w))   # for w=1/n this equals the usual h1=n/(n-1)
 
  # unbiased empirical estimator
  # for w=1/n  the following  would simplify to:  r <- 1/(n-1)*crossprod(xs)
  #r0 <- h1 * t(xs) %*% diag(w) %*% xs
  r0 <- h1 * t(xs) %*% sweep(xs, 1, w, "*") # sweep requires less memory
  
  
  # get ensemble shrinkage intensity
  z <- pvt.get.lambda(xs, lambda, w, protect, verbose=verbose, type="correlation", 0)
  rm(w)
  
  # shrinkage estimate of correlation matrix
  if (z$lambda == 1)
  {
    p <- ncol(xs)
    r <- diag(p)
  }
  else
  {   
    # shrink off-diagonal elements
    r <- (1-z$lambda)*r0
    diag(r) <- 1 
  } 
  rm(xs)
  
  if (z$lambda > 0 && protect > 0) 
  {
    if (verbose)
    {     
      cat("Risk protection: limited translation of the", protect*100, "percent most shrunken components.\n")
    }

    
    # limit the translation/risk of some components
    
    # note: if protect==1 then lambda will be set to zero
    #       if protect==0 full shrinkage will occur
    
    diff <- sm2vec( r0-r )
    rm(r0)
    adiff <- abs(diff)
    sdiff <- sign(diff)
    rm(diff)
  
    M <- quantile(adiff, probs=c(1-protect))
   
    d <- adiff-M # soft thresholding
    rm(adiff)
    d[d < 0] <- 0
    d <- d*sdiff
    rm(sdiff)
  
    W <- vec2sm(d) # correction matrix
    rm(d)
    diag(W) <- 0
    
    r <- r + W  # add correction term
    rm(W)
    
    attr(r, "protect") <- protect
  }
   
  attr(r, "lambda") <- z$lambda
  attr(r, "lambda.estimated") <- z$lambda.estimated
  

  attr(r, "class") <- "shrinkage"

  return(r)   
}    


# compute the *inverse* of the shrinkage correlation matrix

# this procedure exploits the Woodbury identity to efficiently
# compute the inverse of the correlation shrinkage estimator
# directly 

pvt.invscor <- function(x, lambda, w, protect, verbose)
{
  
  if (protect > 0)
  {
    # for protect > 0 the woodbury trick is not applied

    r <- pvt.scor(x=x, lambda=lambda, w=w, protect=protect, verbose=verbose)
    
    ir <- pseudoinverse(r)
    
    attr(ir, "lambda") <- attr(r, "lambda")
    attr(ir, "lambda.estimated") <- attr(r, "lambda.estimated")
    attr(ir, "protect") <- attr(r, "protect")
    rm(r)
    attr(ir, "class") <- "shrinkage"

    return( ir )
  }
  
  
  # standardize input matrix by standard deviations
  xs <- wt.scale(x, w, center=TRUE, scale=TRUE, scale.by="sd") # standardize data matrix

  z <- pvt.get.lambda(xs, lambda, w, protect, verbose=verbose, type="correlation", 0)
  p <- ncol(xs)
    
  if (z$lambda == 1)
  {
    invr <- diag(p)
  }
  else
  {
    # number of zero-variance variables
    zeros <- (attr(xs, "scaled:scale")==0.0)
    
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

  attr(invr, "class") <- "shrinkage"
  
  rownames(invr) <- colnames(xs)
  colnames(invr) <- colnames(xs)
  
  return( invr )
}



# returns lambda to be used for shrinkage
pvt.get.lambda <- function(x, lambda, w, protect, verbose, type=c("correlation", "variance"), target)
{
  type <- match.arg(type)
  
  if (type == "correlation")
  {
     kind <- "lambda (correlation matrix):"
     func <- "C_corlambda"
     
     # note: x needs to be the *scaled* data matrix
  }
  
  if (type == "variance")
  {
     kind <- "lambda.var (variance vector):"
     func <- "C_varlambda"
     
     # note: x needs to be the *centered* data matrix
  }
  
  if (protect == 1) # if all components are translation protected, don't shrink
  {
    lambda <- 0
     cat("Risk protection:  no shrinkage (protect parameter equals 1).\n")
  }
   
  # if lambda = 0: don't shrink
  # if lambda > 0: shrink with given lambda
  # if lambda < 0: shrink with estimated lambda  (the default)
 
     
  if (lambda < 0)
  { 
    if (verbose)
    {
      cat(paste("Estimating optimal shrinkage intensity", kind))     
    }
    
    # estimate optimal shrinkage intensity 
    # target: correlations/covariances -> 0  
    lambda <- .C(func,
            as.double(x),
	    as.integer( nrow(x) ),
	    as.integer( ncol(x) ),
	    as.double(w),
	    as.double(target),
	    lambda=double(1), PACKAGE="corpcor", DUP=FALSE)$lambda

    lambda.estimated <- TRUE
      
    if (verbose)
    {
      cat(paste(" ", round(lambda, 4), "\n", sep=""))     
    }
  }
  else
  {
    if (lambda > 1) lambda <- 1
    lambda.estimated <- FALSE
    
    if (verbose)
    {
      cat(paste("Specified shrinkage intensity", kind, round(lambda, 4), "\n"))     
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




