### cor.shrink.R  (2006-03-09)
###
###    Shrinkage Estimation of Covariance and Correlation Matrix
###    and their inverses
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
   if (missing(lambda)) lambda <- -1  # estimate shrinkage parameter
   w <- pvt.check.w(w, n)
   
   xs <- weighted.scale(x, w)
   r  <- pvt.scor(xs, lambda, w, verbose)

   return(r)
}


# inverse correlation
invcor.shrink <- function(x, lambda, w, verbose=TRUE)
{
   x <- as.matrix(x)
   n <- dim(x)[1]  
   if (missing(lambda)) lambda <- -1  # estimate shrinkage parameter
   w <- pvt.check.w(w, n)
   
   xs <- weighted.scale(x, w)
   invr <- pvt.invscor(xs, lambda, w, verbose)

   return(invr)
}



# covariance
cov.shrink <- function(x, lambda, w, verbose=TRUE)
{   
   x <- as.matrix(x)
   n <- dim(x)[1]   
   if (missing(lambda)) lambda <- -1  # estimate shrinkage parameter
   w <- pvt.check.w(w, n)
   
   xs <- weighted.scale(x, w)
   sc <- attr(xs, "scaled:scale")   
   c <- pvt.scor(xs, lambda, w, verbose)
   c <- sweep(sweep(c, 1, sc, "*"), 2, sc, "*")
                 
   return(c)
}


# precision matrix (inverse covariance)
invcov.shrink <- function(x, lambda, w, verbose=TRUE)
{   
   x <- as.matrix(x)
   n <- dim(x)[1] 
   if (missing(lambda)) lambda <- -1  # estimate shrinkage parameter
   w <- pvt.check.w(w, n)
      
   xs <- weighted.scale(x, w)
   sc <- attr(xs, "scaled:scale")   
   invc <- pvt.invscor(xs, lambda, w, verbose)
   invc <- sweep(sweep(invc, 1, 1/sc, "*"), 2, 1/sc, "*")
   
   return(invc)
}




