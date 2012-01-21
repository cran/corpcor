### shrink.intensity.R  (2012-01-21)
###
###   Functions for computing the shrinkage intensity
###    
###
### Copyright 2005-2012 Julian Schaefer, Rainer Opgen-Rhein and Korbinian Strimmer
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


# estimate shrinkage intensity lambda.var (variance vector)
#
# input:  data matrix
#         weights of each data point
#
estimate.lambda.var = function(x, w, verbose=TRUE)
{
  n = nrow(x) 
  w = pvt.check.w(w, n)
  
  # bias correction factors
  w2 = sum(w*w)           # for w=1/n this equals 1/n   where n=dim(xc)[1]
  h1 = 1/(1-w2)           # for w=1/n this equals the usual h1=n/(n-1)
  h1w2 = w2/(1-w2)        # for w=1/n this equals 1/(n-1)
   
  # center input matrix
  xc = wt.scale(x, w, center=TRUE, scale=FALSE) 

  # compute empirical variances 
  #v = wt.moments(x, w)$var
  v = h1*(colSums(w*xc^2))

  # compute shrinkage target
  target = median(v)


  if (verbose) cat("Estimating optimal shrinkage intensity lambda.var (variance vector): ")

  zz = xc^2
  q1 = colSums( sweep(zz, MARGIN=1, STATS=w, FUN="*") )
  q2 = colSums( sweep(zz^2, MARGIN=1, STATS=w, FUN="*") ) - q1^2   
  numerator = sum( q2 )
  denominator = sum( (q1 - target/h1)^2 )

  if(denominator == 0) 
    lambda.var = 1
  else
    lambda.var = min(1, numerator/denominator * h1w2)
 
  if (verbose) cat(paste(round(lambda.var, 4), "\n")) 
  
  return (lambda.var)
}


# estimate shrinkage intensity lambda (correlation matrix)
#
# input:  data matrix
#         weights of each data point
#
#
# note: the fast algorithm in this function is due to Miika Ahdesm\"aki
#
estimate.lambda = function(x, w, verbose=TRUE)
{
  n = nrow(x)  
  w = pvt.check.w(w, n)
  xs = wt.scale(x, w, center=TRUE, scale=TRUE) # standardize data matrix

  if (verbose) cat("Estimating optimal shrinkage intensity lambda (correlation matrix): ")

  # bias correction factors
  w2 = sum(w*w)           # for w=1/n this equals 1/n   where n=dim(xs)[1]
  h1w2 = w2/(1-w2)        # for w=1/n this equals 1/(n-1)

  sw = sqrt(w)
  Q1.squared = (crossprod(sweep(xs, MARGIN=1, STATS=sw, FUN="*")))^2
  Q2 = crossprod(sweep(xs^2, MARGIN=1, STATS=sw, FUN="*")) - Q1.squared
  denominator = sum(Q1.squared)-sum(diag(Q1.squared)) 
  numerator = sum(Q2)-sum(diag(Q2))

  if(denominator == 0) 
    lambda = 1
  else
    lambda = min(1, numerator/denominator * h1w2)
  
  if (verbose) cat(paste(round(lambda, 4), "\n"))

  return (lambda)
}

