### partial.R  (2006-03-27)
###
###    Partial Correlation and Partial Covariance
###    
###
### Copyright 2003-06 Juliane Schaefer and Korbinian Strimmer
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
cor2pcor <- function(m, tol)
{   
  # invert, then negate off-diagonal entries
  m <- -pseudoinverse(m, tol=tol)
  diag(m) <- -diag(m)

  # standardize and return  
  return(cov2cor(m))
}


#
# backtransformation to correlation matrix
#
# input: partial correlation matrix
# ouput: correlation matrix
pcor2cor <- function(m, tol)
{
  # negate off-diagonal entries, then invert
  m <- -m
  diag(m) <- -diag(m)
  m <- pseudoinverse(m, tol=tol)
  
  # standardize and return 
  return(cov2cor(m))
}



########################################################


cov2pcov <- function(m, tol)
{
  m <- -pseudoinverse(m, tol=tol)
  diag(m) <- -diag(m)

  return(m)
}

pcov2cov <- function(m, tol)
{
  # negate off-diagonal entries, then invert
  m <- -m
  diag(m) <- -diag(m)
  m <- pseudoinverse(m, tol=tol)
  
  return(m)
}


pcor.shrink <- function(x, lambda, w, verbose=TRUE)
{
  pc <- -invcor.shrink(x, lambda, w, verbose)
  diag(pc) <- -diag(pc)

  # standardize and return  
  return(cov2cor(pc))
}


pcov.shrink <- function(x, lambda, lambda.var, w, verbose=TRUE)
{
  pc <- -invcov.shrink(x, lambda, lambda.var, w, verbose)
  diag(pc) <- -diag(pc)

  return(pc)
}


