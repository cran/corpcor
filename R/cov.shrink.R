### cov.shrink.R  (2005-08-10)
###
###    Shrinkage Estimation of Covariance and Correlation Matrix
###
### Copyright 2005 Juliane Schaefer and Korbinian Strimmer
###
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


# compute shrinkage covariance matrix
# 
# 
# S* = lambda T + (1-lambda) S  
#
# S is the unbiased unstructured empirical covariance matrix
# the target T is the unbiased diagonal empirical covariance matrix
#
# if lambda is not specified then it is chosen to minimize the MSE

cov.shrink <- function(x, lambda, verbose=TRUE)
{
  
  vc <- varcov(x, type="unbiased", verbose)
  p <- dim(vc$S)[1]
  
  if (p == 1) return( vc$S )

  if (verbose && p > 50)
    cat("Computing shrinkage covariance matrix\n")
     
     
  # find optimal lambda (for p > 1)
  if (missing(lambda))
  {  
    #offdiagsum.sij.2 <- 0
    #offdiagsum.v.sij <- 0
    #for (i in 1:(p-1))
    #{
    #  for (j in (i+1):p)
    #  {
    #    offdiagsum.sij.2 <- offdiagsum.sij.2 + vc$S[i,j]*vc$S[i,j]
    #    offdiagsum.v.sij <- offdiagsum.v.sij + vc$var.S[i,j]
    #  }
    #}  
       
    offdiagsum.sij.2 <- sum(vc$S[lower.tri(vc$S)]^2)
    offdiagsum.v.sij <- sum(vc$var.S[lower.tri(vc$var.S)])
        
    lambda <- offdiagsum.v.sij/offdiagsum.sij.2
    if (verbose) cat(paste("Estimated shrinkage intensity lambda: ", round(lambda,4), "\n"))
  }
  
  
  #########
  
  if (lambda > 1)
  {
    warning(paste("Overshrinking: intensity lambda set to 1 (allowed range: 0-1)"))
    lambda <- 1  
  }
  if (lambda < 0)
  {
     warning(paste("Undershrinking: intensity lambda set to 0 (allowed range: 0-1)"))
     lambda <- 0  
  }

 
  # construct shrinkage estimator
  S.star <- (1-lambda)*vc$S
  diag(S.star) <- diag(vc$S)
  
  
  attr(S.star, "lambda") <- lambda
  
  
  return( S.star )
}




############### correlations ################




# standardized cov.shrink 
cor.shrink <- function(x, lambda, verbose=TRUE)
{
  return( cov2cor(cov.shrink(x, lambda, verbose)) )
}


# cor2pcor applied to cov.shrink 
pcor.shrink <- function(x, lambda, verbose=TRUE)
{
  return(
    cor2pcor( cov.shrink(x, lambda, verbose),
      exact.inversion=TRUE, check.eigenvalues=FALSE) )
}



