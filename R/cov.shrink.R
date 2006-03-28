### cor.shrink.R  (2006-03-27)
###
###    Shrinkage Estimation of Variance Vector, Correlation Matrix,
###    and Covariance Matrix, and their inverses
###
### Copyright 2005-06 Juliane Schaefer and Korbinian Strimmer
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


# correlation
cor.shrink <- function(x, lambda, w, verbose=TRUE)
{
   x <- as.matrix(x)
   n <- dim(x)[1]
   if (missing(lambda)) lambda <- -1  # estimate correlation shrinkage parameter
   w <- pvt.check.w(w, n)
   
   # shrinkage correlation
   x <- weighted.scale(x, w)
   r <- pvt.scor(x, lambda, w, verbose)

   return(r)
}


# inverse correlation
invcor.shrink <- function(x, lambda, w, verbose=TRUE)
{
   x <- as.matrix(x)
   n <- dim(x)[1]  
   if (missing(lambda)) lambda <- -1  # estimate correlation shrinkage parameter
   w <- pvt.check.w(w, n)
   
   # inverse shrinkage correlation
   wm <- weighted.moments(x, w)
   x <- weighted.scale(x, w, wm=wm)
   invr <- pvt.invscor(wm, x, lambda, w, verbose)

   return(invr)
}


# variances
var.shrink <- function(x, lambda.var, w, verbose=TRUE)
{
   x <- as.matrix(x)
   n <- dim(x)[1]  
   if (missing(lambda.var)) lambda.var <- -1  # estimate variance shrinkage parameter
   w <- pvt.check.w(w, n)
   
   # shrinkage variance
   wm <- weighted.moments(x, w)
   x  <- weighted.scale(x, w, center=TRUE, scale=FALSE, wm=wm) # center data matrix
   sv <- pvt.svar(wm, x, lambda.var, w, verbose)
   
   return(sv)
}


# covariance
cov.shrink <- function(x, lambda, lambda.var, w, verbose=TRUE)
{   
   x <- as.matrix(x)
   n <- dim(x)[1]   
   if (missing(lambda)) lambda <- -1          # estimate correlation shrinkage parameter
   if (missing(lambda.var)) lambda.var <- -1  # estimate variance shrinkage parameter  
   w <- pvt.check.w(w, n)
   
   # shrinkage scale factors
   wm <- weighted.moments(x, w)
   x  <- weighted.scale(x, w, center=TRUE, scale=FALSE, wm=wm) # center data matrix
   sc <- sqrt( pvt.svar(wm, x, lambda.var, w, verbose) )

   # shrinkage correlation
   x <- weighted.scale(x, w, center=FALSE, scale=TRUE, wm=wm) # then standardize data matrix
   c <- pvt.scor(x, lambda, w, verbose)
   
   # shrinkage covariance 
   c <- sweep(sweep(c, 1, sc, "*"), 2, sc, "*")
   attr(c, "lambda.var") <- attr(sc, "lambda.var")
   attr(c, "lambda.var.estimated") <- attr(sc, "lambda.var.estimated")
                    
   return(c)
}


# precision matrix (inverse covariance)
invcov.shrink <- function(x, lambda, lambda.var, w, verbose=TRUE)
{   
   x <- as.matrix(x)
   n <- dim(x)[1] 
   if (missing(lambda)) lambda <- -1          # estimate correlation shrinkage parameter
   if (missing(lambda.var)) lambda.var <- -1  # estimate variance shrinkage parameter  
   w <- pvt.check.w(w, n)

   # shrinkage scale factors
   wm <- weighted.moments(x, w)
   x  <- weighted.scale(x, w, center=TRUE, scale=FALSE, wm=wm) # center data matrix
   sc <- sqrt( pvt.svar(wm, x, lambda.var, w, verbose) )
        
   # inverse shrinkage correlation
   x <- weighted.scale(x, w, center=FALSE, scale=TRUE, wm=wm) # then standardize data matrix
   invc <- pvt.invscor(wm, x, lambda, w, verbose)
   
   # inverse shrinkage covariance 
   invc <- sweep(sweep(invc, 1, 1/sc, "*"), 2, 1/sc, "*")
   attr(invc, "lambda.var") <- attr(sc, "lambda.var")
   attr(invc, "lambda.var.estimated") <- attr(sc, "lambda.var.estimated")
   
   return(invc)
}




