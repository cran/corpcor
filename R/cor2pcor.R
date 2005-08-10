### cor2pcor.R  (2004-09-25)
###
###    Partial Correlation computed by Inversion 
###    of the Covariance or Correlation Matrix
###    
###
### Copyright 2003-04 Juliane Schaefer and Korbinian Strimmer
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


#
# partial correlation matrix
#
# input: covariance matrix or correlation matrix
# ouput: partial correlation matrix
#
cor2pcor <- function(m, exact.inversion=TRUE, check.eigenvalues=TRUE, tol)
{
  if (check.eigenvalues)
  {
    if ( !is.positive.definite(m, tol=tol) )
    {
      stop("Input matrix is not positive definite!")
    }
  }
  
  # standardize
  # m <- cov2cor(m)
  
  # invert, then negate off-diagonal entries
  if (exact.inversion)
  {
    m <- -solve(m)
  }
  else
  {
    m <- -pseudoinverse(m, tol=tol)
  }
  diag(m) <- -diag(m)

  # standardize and return  
  return(cov2cor(m))
}


#
# backtransformation to correlation matrix
#
# input: partial correlation matrix
# ouput: correlation matrix
pcor2cor <- function(m, exact.inversion=TRUE, check.eigenvalues=TRUE, tol)
{
  if (check.eigenvalues)
  {
    if ( !is.positive.definite(m, tol=tol) )
    {
      stop("Input matrix is not positive definite!")
    }
  }

  # standardize
  # m <- cov2cor(m)

  # negate off-diagonal entries, then invert
  m <- -m
  diag(m) <- -diag(m)
  if (exact.inversion)
  {
    m <- solve(m)
  }
  else
  {
    m <- pseudoinverse(m, tol=tol)
  }
  
  # standardize and return 
  return(cov2cor(m))
}


