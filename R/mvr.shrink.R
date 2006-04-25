### mvr.shrink.R  (2005-04-25)
###
###    Fit multivariate linear regression model by shrinkage
###
### Copyright 2006 Rainer Opgen-Rhein and Korbinian Strimmer
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


# compute matrix of regression cofficients
mvr.shrink <- function(x, y, lambda, lambda.var, w, verbose=TRUE)
{
  x <- as.matrix(x)
  n <- nrow(x) # sample size
  p <- ncol(x) # number of predictors
  xnames <- colnames(x)
  if (is.null(xnames)) xnames <- paste("x", 1:p, sep="")

  if (missing(y))
  {
    SEM <- TRUE
    y <- NULL
  }
  else
  {
    SEM <- FALSE
    y <- as.matrix(y)
    
    m <- ncol(y) # dimension of response
    if (nrow(y) != n) stop("Sample size must be identical in x and y.")
    
    ynames <- colnames(y)
    if (is.null(ynames)) ynames <- paste("y", 1:m, sep="")
  }
  
  
  # mean and partical covariance of combined data
  yx <- cbind(y, x) 
  mu <- wt.moments(yx, w)$mean
  pcov <- pcov.shrink(yx, lambda=lambda, lambda.var=lambda.var, w=w, verbose=verbose) 
  
  # estimate regression coefficients
  if (SEM) # consider each predictor variable in turn as response 
  {
    beta <- matrix(NA, ncol=p+1, nrow=p)
    for (i in 1:p)
    {
      beta[i,-i-1] <- pvt.usr(i, -i, pcov, mu)
      beta[i, i+1] <- 0
    }  
    rownames(beta) <- xnames
  }
  else # use specified response variables
  {
    beta <- matrix(NA, ncol=p+1, nrow=m)
    for (i in 1:m)
      beta[i,] <- pvt.usr(i, -(1:m), pcov, mu)
    rownames(beta) <- ynames
  }
  
  colnames(beta) <- c("(Intercept)", xnames)
  attr(beta,"lambda") <- attr(pcov,"lambda")
  attr(beta,"lambda.estimated") <- attr(pcov,"lambda.estimated")
  attr(beta,"lambda.var") <- attr(pcov,"lambda.var")
  attr(beta,"lambda.var.estimated") <- attr(pcov,"lambda.var.estimated")
  
  attr(beta,"class") <- "shrinkage"

  return( beta )
}


# predict response
mvr.predict <- function(coef, x)
{
  y <- as.matrix(x) %*% t(coef[,-1, drop=FALSE])  
  y <- sweep(y, 2, coef[,1, drop=FALSE], "+")

  rownames(y) <- rownames(x)

  return(y)
}


##### private function #######

# unvariate shrinkage regression
pvt.usr <- function(response, predictors, pcov, mu)
{
   b <- pcov[response,predictors]/pcov[response,response]
   a <- mu[response] - sum(b * mu[predictors]) # intercept
   
   return ( c(a, b) )
}

