### shrink.internal.R  (2008-12-01)
###
###    Non-public functions used in the covariance shrinkage estimator 
###    
###
### Copyright 2005-08 Korbinian Strimmer
###
### This file is part of the `corpcor' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
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
pvt.check.w = function(w, n)
{
   if (missing(w))
   {
     w = rep(1/n, n)   # equal weights
   }
   else
   {
     if (length(w) != n)
     {
       stop("Weight vector has incompatible length", call. = FALSE)
     } 
     
     sw = sum(w)
     if (sw != 1) 
       w = w/sw      # make weights sum up to 1
   } 
 
   return(w)
}


# print function
print.shrinkage = function(x, ...)
{
  attr(x, "class") = NULL
  
  lambda = attr(x, "lambda")
  lambda.estimated = attr(x, "lambda.estimated")
  attr(x, "lambda") = NULL 
  attr(x, "lambda.estimated") = NULL 
  
  lambda.var = attr(x, "lambda.var")
  lambda.var.estimated = attr(x, "lambda.var.estimated")
  attr(x, "lambda.var") = NULL 
  attr(x, "lambda.var.estimated") = NULL 
  
  spv = attr(x, "spv")
  attr(x, "spv") = NULL

  NextMethod("print", x, quote = FALSE, right = TRUE)
    
  cat("\n")
  if (!is.null(lambda.estimated))  
  {
    if (lambda.estimated) 
      le = "(estimated)"
    else
      le = "(specified)"  
    cat(paste("Shrinkage intensity lambda (correlation matrix):", round(lambda,4), le, "\n"))
  }
  
  if (!is.null(lambda.var.estimated))  
  {
    if (lambda.var.estimated)
      lve = "(estimated)"
    else
      lve = "(specified)"  
    cat(paste("Shrinkage intensity lambda.var (variance vector):", round(lambda.var, 4), lve, "\n"))
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

pvt.svar = function(x, lambda.var, w, verbose)
{    
  # center input matrix
  xs = wt.scale(x, w, center=TRUE, scale=FALSE) 
  
  # compute variances 
  v = wt.moments(xs, w)$var
  
  # compute target
  tgt = median(v)
          
  # shrinkage estimate 
  z = pvt.get.lambda(xs, lambda.var, w, verbose=verbose, type="variance", tgt)   
  vs = z$lambda.var*tgt + (1-z$lambda.var)*v
         
  attr(vs, "lambda.var") = z$lambda.var
  attr(vs, "lambda.var.estimated") = z$lambda.var.estimated
  attr(vs, "class") = "shrinkage"
  
  return(vs)   
}    


# returns lambda to be used for shrinkage
pvt.get.lambda = function(x, lambda, w, verbose, type=c("correlation", "variance"), target)
{
  type = match.arg(type)
  
  if (type == "correlation")
  {
     kind = "lambda (correlation matrix):"
     func = "C_corlambda"
     
     # note: x needs to be the *scaled* data matrix
  }
  
  if (type == "variance")
  {
     kind = "lambda.var (variance vector):"
     func = "C_varlambda"
     
     # note: x needs to be the *centered* data matrix
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
    lambda = .C(func,
            as.double(x),
	    as.integer( nrow(x) ),
	    as.integer( ncol(x) ),
	    as.double(w),
	    as.double(target),
	    lambda=double(1), PACKAGE="corpcor", DUP=FALSE)$lambda

    lambda.estimated = TRUE
      
    if (verbose)
    {
      cat(paste(" ", round(lambda, 4), "\n", sep=""))     
    }
  }
  else
  {
    if (lambda > 1) lambda = 1
    lambda.estimated = FALSE
    
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


