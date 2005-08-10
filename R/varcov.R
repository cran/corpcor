### cov.shrink.R  (2005-06-07)
###
###    Variance of the Covariance  Matrix
###
### Copyright 2005 Korbinian Strimmer
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


# compute the empirical covariance matrix S=cov(x) given a data matrix x
# as well as the *variances* associated with the individual entries S[i,j]
#

varcov <- function(x, type=c("unbiased", "ML"), verbose=TRUE)
{
    if (!is.matrix(x)) x <- as.matrix(x)     
    n <- dim(x)[1]
    p <- dim(x)[2]
 
    if (verbose && p > 50)
      cat("Computing empirical covariance matrix and its variance\n")
       
     
    # weights for the "unbiased" and "ML" cases
    type <- match.arg(type)
    if (type=="unbiased")
    {
      h1 <- 1/(n-1)
      h2 <- n/(n-1)/(n-1)
    }    
    if (type=="ML")
    {
      h1 <- 1/n
      h2 <- (n-1)/n/n
    }
 
    s <- matrix(NA, ncol=p, nrow=p)   
    vs <- matrix(NA, ncol=p, nrow=p)
    xc <- scale(x, scale=FALSE) # center the data
    
    # diagonal elements
    for (i in 1:p)
    {
      zii <- xc[,i]^2
      s[i,i] <- sum(zii)*h1
      vs[i,i] <- var(zii)*h2
    }
    
    if (p == 1) return(list(S=s, var.S=vs))
    
    if (verbose && p > 50)
      cat(paste("Wait for", p, "points (50 per row):\n")) 
    
    # off-diagonal elements
    for (i in 1:(p-1))
    {
      if (verbose && p > 50)
      {
        cat(".")
        if (i %% 50 == 0) cat(paste(" ", i, "\n"))
      }
      
      for (j in (i+1):p)
      {
        zij <- xc[,i]*xc[,j] 
	s[i,j] <- sum(zij)*h1
        s[j,i] <- s[i,j]
        
        vs[i,j] <- var(zij)*h2
        vs[j,i] <- vs[i,j]	 
      }
      
    }
    if (verbose && p > 50) cat(paste(". ", i+1, "\n"))

    return(list(S=s, var.S=vs))
}

