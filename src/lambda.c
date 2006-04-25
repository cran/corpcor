
/* lambda.c  (2006-04-24)  
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
   
  /* variance of mean */
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
  double h2, sw2;
  
  /* bias correction factors */
  sw2 = 0;
  for (i=0; i < n; i++)
    sw2 += w[i]*w[i];
  
  *h1 = 1/(1-sw2);           /* for w=1/n this equals the usual h1=n/(n-1) */
   h2 = (*h1) * (*h1)*sw2;   /* for w=1/n this equals h2=n/(n-1)^2         */
  *h3 = (*h1) * h2;          /* for w=1/n this equals h3=n^2/(n-1)^3       */ 
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
 * input:  lr      limit risk of individual components
 *         xs      scaled matrix 
 *         n       sample size (number of rows)
 *         p       number of variables (number of columns)
 *         w       data weights
 * output: lambda  shrinkage intensity (Target: correlations -> 0)  
 */
void C_corlambda(int *lr, double* xs, int* n, int* p, double* w, double* lambda)
{
  double h1, h3, rr, vr;
  double a, b, suma, sumb;
  int i, k, l, nn, pp, kn, ln; 
  double* xsij;
  double la, minla;
  
  nn = *n;
  pp = *p;

  /* allocate vector - error handling is done by R */
  xsij = (double *) Calloc((size_t) nn, double);

 
  /* compute numerator and denominator for the optimal 
     shrinkage intensity using target: zero off-diagonal entries */

  la = 0;     /* component-wise lambda */
  minla = 1;  /* smallest component-wise lambda */
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

       C_compute_lambda(a, b, &la);
       if (la < minla) /* find smallest component lambda */
         minla = la;
    }
  }
  
 
  /* optimal ensemble shrinkage intensity */
  C_compute_lambda(suma, sumb, lambda);
  
  /* make sure that MSE of every individual component always decreases */
  if (*lr == 1) /* limit.risk==TRUE */
    if ( *lambda > 2*minla) *lambda = 2*minla;

  
  /* free vector */
  Free(xsij);  
}




/* 
 * input:  lr      limit risk of individual components
 *         xc      centered matrix 
 *         n       sample size (number of rows)
 *         p       number of variables (number of columns)
 *         w       data weights
 * output: lambda  shrinkage intensity (Target: average empirical variance)  
 */
void C_varlambda(int* lr, double* xc, int* n, int* p, double* w, double* lambda)
{
  double h1, h3, vv, varvv, fac;
  double suma, sumb,meanb, la, minla;
  int i, k, nn, pp, kn; 
  double* xc2, *avec, *bvec;
    
  nn = *n;
  pp = *p;

  C_biascorrectionfactor(w, nn, &h1, &h3);


  /* allocate vectors - error handling is done by R */
  xc2 = (double *) Calloc((size_t) nn, double);
  avec = (double *) Calloc((size_t) pp, double);
  bvec = (double *) Calloc((size_t) pp, double);
 
    
  /* compute numerator and denominator for the optimal 
     shrinkage intensity using target: mean of the varianes */
 
  vv = 0;
  varvv = 0;
  meanb = 0;
  fac = (1.0-1.0/((double)pp));
  for (k=0; k < pp; k++)
  {  
       kn = k*nn;
  
       /* xc[,k] * xc[,k] */
       for (i=0; i < nn; i++)
       {
           xc2[i] = xc[kn+i] * xc[kn+i];       
       }
 
       C_meanvarmean(xc2, nn, w,  &vv, &varvv);
              
       avec[k] = fac*h3*varvv;
       bvec[k] = h1*vv;
       meanb += bvec[k];
  }
  meanb = meanb/pp;
  
  la = 0;     /* component-wise lambda */
  minla = 1;  /* smallest component-wise lambda */
  suma = 0;
  sumb = 0;
  for (k=0; k < pp; k++)
  {  
    bvec[k] = bvec[k] - meanb;
    bvec[k] = bvec[k]*bvec[k];
  
    suma += avec[k];
    sumb += bvec[k];
    
    C_compute_lambda(avec[k], bvec[k], &la);
    if (la < minla) /* find smallest component lambda */
      minla = la;
  }
  
  /* optimal ensemble shrinkage intensity */
  C_compute_lambda(suma, sumb, lambda);
  
  /* make sure that MSE of every individual component always decreases */
  if (*lr == 1) /* limit.risk==TRUE */
    if ( *lambda > 2*minla) *lambda = 2*minla;
  
     
  /* free vectors */
  Free(xc2);
  Free(avec);
  Free(bvec);
}
