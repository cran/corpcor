
/* lambda.c  (2006-07-13)  
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
          

/* estimate weighted mean and weighted variance from data vector z*/

void C_meanvarmean(double* z, int n, double* w, double* m, double* v)
{
  int i;
  double zc, mean, varmean;
  
  /* weighted mean */
  mean = 0;
  for (i=0; i < n; i++)
  {
    mean += w[i]*z[i];
  }  
   
  /* weighted variance */
  varmean = 0;
  for (i=0; i < n; i++)
  {
    zc = z[i]-mean;
    varmean += w[i]*zc*zc;
  }
  
  /* return */
  *m = mean;
  *v = varmean; 
}


void C_biascorrectionfactor(double *w, int n, double* h1, double* h3)
{
  int i;
  double sw2;
  
  /* bias correction factors */
  sw2 = 0;
  for (i=0; i < n; i++)
    sw2 += w[i]*w[i];
  
  *h1 = 1/(1-sw2);              /* for w=1/n this equals  n/(n-1)        */ 
  *h3 = (*h1)*(*h1)*(*h1) *sw2; /* for w=1/n this equals  n^2/(n-1)^3    */ 
}



void C_compute_lambda(double numerator, double denominator, double* lambda) 
{

  if (denominator == 0.0)
  {
    /* in case of equality of target and unconstrained estimate
       (= zero misspecification) set shrinkage intensity to 1 */
    
    *lambda = 1;
  }
  else
  {
    *lambda = numerator/denominator;
    
    /* don't overshrink */
    if( (*lambda) > 1 ) *lambda = 1;
  }
}



/* 
 * input:  xs      scaled matrix 
 *         n       sample size (number of rows)
 *         p       number of variables (number of columns)
 *         w       data weights
 *         t       target (ignored: in this function target is always assumed to be zero)      
 * output: lambda  shrinkage intensity  
 */
void C_corlambda(double* xs, int* n, int* p, double* w, double* t, double* lambda)
{
  double h1, h3, rr, vr;
  double a, b, suma, sumb;
  int i, k, l, nn, pp, kn, ln; 
  double* xsij;
  
  nn = *n;
  pp = *p;

  /* allocate vector - error handling is done by R */
  xsij = (double *) Calloc((size_t) nn, double);

 
  /* compute numerator and denominator for the optimal 
     shrinkage intensity using target: zero off-diagonal entries */

  suma = 0;
  sumb = 0;
  rr = 0;
  vr = 0;

  C_biascorrectionfactor(w, nn, &h1, &h3);
 
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
 
       C_meanvarmean(xsij, nn, w, &rr, &vr);
        
       a = vr*h3;
       b = (rr*h1);
       b = b*b;
       suma += a;
       sumb += b;
    }
  }
  
 
  /* optimal ensemble shrinkage intensity */
  C_compute_lambda(suma, sumb, lambda);
  
  
  /* free vector */
  Free(xsij);  
}




/* 
 * input:  xc      centered matrix 
 *         n       sample size (number of rows)
 *         p       number of variables (number of columns)
 *         w       data weights
 *         t       target (this is the median value of all empirical variances)      
 * output: lambda  shrinkage intensity (Target: average empirical variance)  
 */
void C_varlambda(double* xc, int* n, int* p, double* w, double* t, double* lambda)
{
  double h1, h3, vv, varvv;
  double a, b, suma, sumb;
  int i, k, nn, pp, kn; 
  double* xc2;
    
  nn = *n;
  pp = *p;

  C_biascorrectionfactor(w, nn, &h1, &h3);


  /* allocate vectors - error handling is done by R */
  xc2 = (double *) Calloc((size_t) nn, double);
 
    
  /* compute numerator and denominator for the optimal 
     shrinkage intensity using the target t */
 
  vv = 0;
  varvv = 0;
  suma = 0;
  sumb = 0;
  
  for (k=0; k < pp; k++)
  {  
       kn = k*nn;
  
       /* xc[,k] * xc[,k] */
       for (i=0; i < nn; i++)
       {
           xc2[i] = xc[kn+i] * xc[kn+i];       
       }
 
       C_meanvarmean(xc2, nn, w,  &vv, &varvv);
              
       a = h3*varvv;
       b = h1*vv - (*t); 
       suma += a;
       sumb += b*b;
  }

  
  /* optimal ensemble shrinkage intensity */
  C_compute_lambda(suma, sumb, lambda);
  
     
  /* free vectors */
  Free(xc2);
}
