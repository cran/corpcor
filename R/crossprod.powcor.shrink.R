### crossprod.powcor.shrink.R  (2011-06-27)
###
###    Efficient computation of crossprod(R^alpha, y)
###
### Copyright 2011 A. Pedro Duarte Silva and Korbinian Strimmer
###
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


# computes R_shrink^alpha %*% y
crossprod.powcor.shrink = function(x, y, alpha, lambda, w, verbose=TRUE)
{
  if (missing(alpha)) stop("Please specify the exponent alpha!")

   x = as.matrix(x)
   y = as.matrix(y)
   p = ncol(x)
   if (nrow(y) != p) stop("Matrix y must have ", p, " rows!")
   n = nrow(x)  
   if (missing(lambda)) lambda = -1  # estimate correlation shrinkage parameter
   w = pvt.check.w(w, n)
   
   # crossprod of matrix power of shrinkage correlation with y
   cp.powr = pvt.cppowscor(x=x, y=y, alpha=alpha, lambda=lambda, w=w, verbose=verbose)
   if (verbose) cat("\n")
   
   return(cp.powr)
}



##### internal functions ######

# this procedure exploits a special identity to efficiently
# compute the crosspord of matrix power of the correlation shrinkage estimator with y

# computes R_shrink^alpha %*% y
pvt.cppowscor = function(x, y, alpha, lambda, w, verbose)
{  
  # standardize input matrix by standard deviations
  xs = wt.scale(x, w, center=TRUE, scale=TRUE) # standardize data matrix

  # bias correction factor
  h1 = 1/(1-sum(w*w))   # for w=1/n this equals the usual h1=n/(n-1)

  z = pvt.get.lambda(xs, lambda, w, verbose=verbose, type="correlation", 0)
  p = ncol(xs)
    
  if (z$lambda == 1 | alpha == 0) # result in both cases R is the identity matrix
  {
      cp.powr = y    # return y
  }
  else
  {
    # number of zero-variance variables
    zeros = (attr(xs, "scaled:scale")==0.0)
    
    svdxs = fast.svd(xs)
    m = length(svdxs$d)  # rank of xs
               
    UTWU = t(svdxs$u) %*% sweep(svdxs$u, 1, w, "*") #  t(U) %*% diag(w) %*% U
    C = sweep(sweep(UTWU, 1, svdxs$d, "*"), 2, svdxs$d, "*") # D %*% UTWU %*% D
    C = (1-z$lambda) * h1 * C

    C = (C + t(C))/2  # symmetrize for numerical reasons (mpower() checks symmetry)
    
    # note: C is of size m x m, and diagonal if w=1/n
         
    if (lambda==0.0) # use eigenvalue decomposition computing the matrix power
    {
      if (m < p-sum(zeros)) 
        warning(paste("Estimated correlation matrix doesn't have full rank",
          "- pseudoinverse used for inversion."), call. = FALSE)    
       
      cp.powr =  svdxs$v %*%  (mpower(C, alpha)  %*% crossprod( svdxs$v, y))
    }
    else # use a special identity for computing the matrix power
    {
      F = diag(m) - mpower(C/z$lambda + diag(m), alpha)
      cp.powr = (y - svdxs$v %*% (F %*% crossprod(svdxs$v, y) ))*(z$lambda)^alpha 
    }
 
    # set all diagonal entries in R_shrink corresponding to zero-variance variables to 1
    cp.powr[zeros,] = y[zeros,]
  }
  rownames(cp.powr) = colnames(xs)
  colnames(cp.powr) = colnames(y)  
  rm(xs)

  attr(cp.powr, "lambda") = z$lambda
  attr(cp.powr, "lambda.estimated") = z$lambda.estimated
  attr(cp.powr, "class") = "shrinkage"
  
  return( cp.powr )
}




