### cor.shrink.R  (2005-09-28)
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



cor.shrink <- function(x, lambda, verbose=TRUE)
{
  sx <- scale(x)  # standardize data (and turn x into a matrix)
 
  p <- dim(sx)[2]
  if (p == 1) return( as.matrix(1) ) 
  
  # estimate variance of empirical correlation coefficients 
  vc <- varcov(sx, type="unbiased", verbose)
  
  # find optimal lambda
  if (missing(lambda))
  {   
    offdiagsum.rij.2 <- sum(vc$S[lower.tri(vc$S)]^2)
    offdiagsum.v.rij <- sum(vc$var.S[lower.tri(vc$var.S)])
        
    lambda <- offdiagsum.v.rij/offdiagsum.rij.2
    if (verbose) cat(paste("Estimated shrinkage intensity lambda: ",
        round(lambda,4), "\n"))
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
  
  #########
 
  # construct shrinkage estimator
  R.star <- (1-lambda)*vc$S
  diag(R.star) <- rep(1, p)
  
  attr(R.star, "lambda") <- lambda
  
  return( R.star )
}


############### derived estimators ################


cov.shrink <- function(x, lambda, verbose=TRUE)
{
   if( !is.matrix(x) ) x <- as.matrix(x)

   # shrinkage correlation coefficients
   R.star <- cor.shrink(x, lambda=lambda, verbose=verbose)

   # unbiased empirical variances
   V <- apply(x, 2, var)
     
   return( rebuild.cov(R.star, V) )
}



# cor2pcor applied to cor.shrink 
pcor.shrink <- function(x, lambda, verbose=TRUE)
{
  return(
    cor2pcor( cor.shrink(x, lambda=lambda, verbose=verbose),
      exact.inversion=TRUE, check.eigenvalues=FALSE) )
}



