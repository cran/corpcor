
/* lambda.c  (2005-03-27)  
 *
 * Copyright 2006 Korbinian Strimmer
 *
 * This file is part of the `corpcor' library for R and related languages.
 * It is made available under the terms of the GNU General Public
 * License, version 2, or at your option, any later version,
 * incorporated herein by reference.
 *
 * This program is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details. 
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA 
 */


#include <R.h>
          

/* estimate weighted mean and variance of mean from data vector z*/

void C_meanvarmean(double* z, int n, double* w, double* h1, double* h3, 
     double* m, double* v)
{
  int i;
  double zc, mean, varmean;
  

  /* weighted mean */
  mean = 0;
  for (i=0; i < n; i++)
  {
    mean += w[i]*z[i];
  }
  mean = (*h1)*mean;
  
   
  /* variance of mean */
  varmean = 0;
  for (i=0; i < n; i++)
  {
    zc = z[i]-mean;
    varmean += w[i]*zc*zc;
  }
  varmean = (*h3)*varmean;
  
  /* return */
  *m = mean;
  *v = varmean; 
}


void C_biascorrectionfactor(double *w, int n, double* h1, double* h3)
{
  int i;
  double h2, sw2;
  
  /* bias correction factors */
  sw2 = 0;
  for (i=0; i < n; i++)
    sw2 += w[i]*w[i];
  
  *h1 = 1/(1-sw2);           /* for w=1/n this equals the usual h1=n/(n-1) */
   h2 = (*h1) * (*h1)*sw2;   /* for w=1/n this equals h2=n/(n-1)^2         */
  *h3 = (*h1) * h2;          /* for w=1/n this equals h3=n^2/(n-1)^3       */ 
}





/* 
 * input:  xs      scaled matrix 
 *         n       sample size (number of rows)
 *         p       number of variables (number of columns)
 *         w       data weights
 * output: lambda  shrinkage intensity (Target: correlations -> 0)  
 */
void C_corlambda(double* xs, int* n, int* p, double* w, double* lambda)
{
  double h1, h3, rr, vr;
  double numerator, denominator;
  int i, k, l, nn, pp, kn, ln; 
  double* xsij;
  
  nn = *n;
  pp = *p;

  C_biascorrectionfactor(w, nn, &h1, &h3);


  /* allocate vector - error handling is done by R */
  xsij = (double *) Calloc((size_t) nn, double);

 
  /* compute numerator and denominator for the optimal 
     shrinkage intensity using target: zero off-diagonal entries */

  denominator = 0;
  numerator = 0;
  rr = 0;
  vr = 0;
  
  for (k=0; k < pp-1; k++)
  {
    for (l=k+1; l < pp; l++)
    {
       kn = k*nn;
       ln = l*nn;
  
       /* xs[,k] * xs[,l] */
       for (i=0; i < nn; i++)
       {
           xsij[i] = xs[kn+i] * xs[ln+i];       
       }
 
       C_meanvarmean(xsij, nn, w, &h1, &h3,  &rr, &vr);
        
       numerator += vr;
       denominator += rr*rr;
    }
  }

  
  /* return estimated shrinkage intensity */
   
  if (denominator == 0.0)
  {
    /* in case of equality of target and unconstrained estimate
       (= zero misspecification) set shrinkage intensity to 1 */
    
    *lambda = 1;
  }
  else
  {
    /* optimal shrinkage intensity for target with zero offdiagonal elements */
    *lambda = numerator/denominator;
    
    /* don't overshrink */
    if( (*lambda) > 1 ) *lambda = 1;
  }
  
  /* free vector */
  Free(xsij);  
}




/* 
 * input:  xc      centered matrix 
 *         n       sample size (number of rows)
 *         p       number of variables (number of columns)
 *         w       data weights
 * output: lambda  shrinkage intensity (Target: average empirical variance)  
 */
void C_varlambda(double* xc, int* n, int* p, double* w, double* lambda)
{
  double h1, h3, vv, varvv, meanvv, tmp;
  double numerator, denominator;
  int i, k, nn, pp, kn; 
  double* xc2;
  double* varvec;
    
  nn = *n;
  pp = *p;

  C_biascorrectionfactor(w, nn, &h1, &h3);


  /* allocate vectors - error handling is done by R */
  xc2 = (double *) Calloc((size_t) nn, double);
  varvec = (double *) Calloc((size_t) pp, double);
 
    
  /* compute numerator and denominator for the optimal 
     shrinkage intensity using target: mean of the varianes */
 
  vv = 0;
  varvv = 0;
  meanvv = 0;
  numerator = 0;
  denominator = 0;

  for (k=0; k < pp; k++)
  {  
       kn = k*nn;
  
       /* xc[,k] * xc[,k] */
       for (i=0; i < nn; i++)
       {
           xc2[i] = xc[kn+i] * xc[kn+i];       
       }
 
       C_meanvarmean(xc2, nn, w, &h1, &h3,  &vv, &varvv);
              
       varvec[k] = vv;
       meanvv += vv;
       numerator += varvv;
  }
  meanvv = meanvv/pp;
    
  numerator = (1.0-1.0/((double)pp))*numerator;
  denominator = 0;   
  for (k=0; k < pp; k++)
  {  
     tmp = varvec[k] - meanvv;
     denominator += tmp*tmp;
  }
      
  /* return estimated shrinkage intensity */
   
  if (denominator == 0.0)
  {
    /* in case of equality of target and unconstrained estimate
       (= zero misspecification) set shrinkage intensity to 1 */
    
    *lambda = 1;
  }
  else
  {
    /* optimal shrinkage intensity */
    *lambda = numerator/denominator;
    
    /* don't overshrink */
    if( (*lambda) > 1 ) *lambda = 1;
  }
  
  
  /* free vectors */
  Free(xc2);
  Free(varvec);
}

