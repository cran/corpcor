### cov.bagged.R  (2004-03-15)
###
###     Variance reduced estimators of cov, cor, and pcor
###     using bootstrap aggregation ("bagging")
###
### Copyright 2003-04 Juliane Schaefer and Korbinian Strimmer
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




# bagged estimators

cov.bagged <- function(x, R=1000, ...)
{
  vec.out <- bag.fun(cov, x, R=R, diag=TRUE, ...)
  mat.out <- vec2sm(vec.out, diag=TRUE)
  
  return( mat.out )
}

cor.bagged <- function(x, R=1000, ...)
{
  vec.out <- bag.fun(cor, x, R=R, diag=FALSE, ...)
  mat.out <- vec2sm(vec.out, diag=FALSE)
  diag(mat.out) <- rep(1, dim(mat.out)[1]) # fill diagonal with 1
  
  return( mat.out )
}

pcor.bagged <- function(x, R=1000, ...)
{
  vec.out <- bag.fun(pcor.pseudo, x, R=R, diag=FALSE, ...)
  mat.out <- vec2sm(vec.out, diag=FALSE)
  diag(mat.out) <- rep(1, dim(mat.out)[1]) # fill diagonal with 1
  
  return( mat.out )
}


#############################################################################

# internal

bag.fun <- function(fun, data, R, diag, ...)
{
  # number of variables 
  p <- dim(data)[2]
  
  # index vector for lower triangle
  lo <- lower.tri(matrix(NA, nrow=p, ncol=p), diag=diag)

  # bootstrap function
  boot.fun <- function(data, i) 
  {
    vec <- as.vector( fun(data[i,], ...)[lo] )
      
    # if we get NAs flag result as being erroneous
    if (sum(is.na(vec)) > 0) class(vec) <- "try-error"

    return( vec )
  }   
     
  #bag variable 
  #boot.out <- boot(data=data, statistic=boot.fun, R=R)
  boot.out <- robust.boot(data=data, statistic=boot.fun, R=R)
  
  bag <- apply( boot.out$t, 2, mean)
    
  return( bag )
}


# simple bootstrap function (robust against errors)
robust.boot <- function(data, statistic, R)
{
  idx <- 1:dim(data)[1]
  
  # determine dimension of statistic
  repeat
  {
    bx <- sample(idx, replace=TRUE)
    val <- try(statistic(data, bx)) 
    
    if (class(val) != "try-error") break
  }
  dim.statistic <- length(val)
  output <- matrix(nrow=R, ncol=dim.statistic)
  
  replicate.count <- 0
  error.count <- 0
  while (replicate.count < R)
  {
    bx <- sample(idx, replace=TRUE)
    val <- try(statistic(data, bx)) 
    
    if (class(val) == "try-error") # if we get a numerical error we simply repeat the draw ..
    {
      error.count <- error.count+1
      #cat("Bootstrapping continues, drawing an alternative bootstrap sample ...\n")
      
      if (error.count > R) stop("Too many errors encountered during the bootstrap.")
    }
    else
    {
      replicate.count <- replicate.count+1
      output[replicate.count,] <- val
    }
  }
  
  if (error.count > 0) warning(paste(error.count, "out of", R,
   "bootstrap samples were repeated due to errors."))
  
  return(list(t=output))
} 



# for pcor.bagged


#
# compute partial correlations given the data x 
# using the pseudoinverse of cor (x)
#

pcor.pseudo <- function(x, tol)
{
  pc <- -psinv.cor(x, tol)
  diag(pc) <- -diag(pc)
  
  return(cov2cor(pc)) 
}


#
# compute inverse correlation matrix
#
# this is numerically equivalent to pseudoinverse(cor(x)) but much
# faster for n << p 
#
psinv.cor <- function (x, tol)
{
    n <- dim(x)[1]
    p <- dim(x)[2]
    
    xs <- scale(x) # standardize
   
    if (n < p)
    {
        # the following is *much* faster than inverting the
	# p x p cor matrix directly
         
        xsvd <- fast.svd(xs, tol)   # fast svd on "fat" matrix (using svd on n x n matrix)
        if (length(xsvd$d) == 0)
        {
           ic <- array(0, c(p,p))
        }
        else
        {
           ic <- xsvd$v %*% ((n-1)/xsvd$d^2 * t(xsvd$v))
        }
        
    }
    else
    {
	ic <- pseudoinverse(crossprod(xs)/(n-1), tol)   # invert p x p matrix using svd	
    }
      
    return(ic)
}

