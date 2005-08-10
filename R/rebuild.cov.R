### rebuild.cov.R (2004-01-15)
###
###    Rebuild Covariance Matrix from Correlation Matrix and Variances
###    
###
### Copyright 2003-04 Korbinian Strimmer
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
# rebuild covariance matrix
#
# input:  correlation matrix           rho(ij)
#         vector with variances        var(i) 
# output: correlation matrix   rho(ij)*sqrt(var(i)*var(j))
#
rebuild.cov <- function(r, v)
{
  resid.sd <- sqrt(v)
  m <- sweep(sweep(r, 1, resid.sd, "*"), 2, resid.sd, "*") 
    
  return(m)
}


