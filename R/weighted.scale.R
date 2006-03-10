### weighted.scale.R  (2006-03-09)
###
###    Weighted Expectations and Variances
###    
###
### Copyright 2006 Rainer Opgen-Rhein and Korbinian Strimmer
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


# in all the following functions, 
# w is a vector of weights with sum(w)=1


# mean and variance


# this exists already in R 
#weighted.mean <- function(xvec, w)
#{
#  return( sum(w*xvec) )
#}

weighted.var <- function(xvec, w) 
{
  w <- pvt.check.w(w, length(xvec))

  # bias correction factor
  h1 <- 1/(1-sum(w*w))   # for w=1/n this equals the usual h1=n/(n-1)
       
  xc <- xvec-weighted.mean(xvec, w)
  s2 <- h1*weighted.mean(xc*xc, w)

  return( s2 ) 
}


weighted.moments <- function(x, w)
{
  x <- as.matrix(x)
  w <- pvt.check.w(w, dim(x)[1])
     
  m <- apply(x, 2, weighted.mean, w=w)
  v <- apply(x, 2, weighted.var, w=w)
  
  return( list(mean=m, var=v) )
}


# scale using the weights
weighted.scale <- function(x, w, center=TRUE, scale=TRUE)
{
  x <- as.matrix(x)
  w <- pvt.check.w(w, dim(x)[1])
  
  # compute column means and variances
  wm <- weighted.moments(x, w)

  if (center==TRUE)
  {
      x <- sweep(x, 2, wm$mean, "-")	
      attr(x, "scaled:center") <- wm$mean
  } 
  
  if (scale==TRUE)
  {
      sd <- sqrt(wm$var)
      x <- sweep(x, 2, sd, "/")
      attr(x, "scaled:scale") <- sd
  } 

  return(x)
}



