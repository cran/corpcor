### mpower.R  (2008-11-30)
###
###    Compute the Power of a Real Symmetric Matrix
###
### Copyright 2008 Korbinian Strimmer
###
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



# compute m^alpha where m is a matrix
mpower = function(m, alpha)
{
  em = eigen(m, symmetric=TRUE)
  ma = em$vectors %*% tcrossprod(diag(em$values^alpha), em$vectors)
 
  return( ma ) 
}
