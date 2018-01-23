/*
 * Copyright 2018 Hatim Chergui, Mustapha Benjillali, 
 * and Mohamed-Slim Alouini. Contact email <chergui@ieee.org>
 * ***************************************************************************************************************
 * Citation: If you use this software or any (modified) part of it, please cite it as:
 * Hatim Chergui, Mustapha Benjillali and Mohamed-Slim Alouini. (2018, January 22). 
 * Multivariate Fox H-Function C/MEX Package: mfoxh (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.1157194
 * **************************************************************************************************************
 *
 * The quasi-Monte Carlo (QMC) complex integration has been devloped by 
 * extending the online real-valued QMC module https://github.com/diazona 
 * to the complex domain.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_matrix_int.h>
#include "mfox.h"
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
    
    
    size_t i, j, k, mx, nx, dim;
    double  *ind, *xr, *xi, *zr, *zi, *max_call, *tol;
    gsl_complex x;
    gsl_vector_complex *xl, *xu;
    gsl_matrix_int *index;
    gsl_matrix_complex *Arg[20], *Emp;

if((nrhs -7) % 2  == 0){
dim = (size_t)((nrhs -7)/2);}
else{mexErrMsgTxt("Number of inputs is incorrect\n");}

// Retrieve vector [0, n, m1, n1, ..., mM, nM]
if(mxIsComplex(prhs[0])){mexErrMsgTxt("Indices should be integers\n");}
nx = mxGetN(prhs[0]);
if(nx != 2*dim+2){mexErrMsgTxt("Missing input(s) or extra elements in the first argument (index)\n");}
index   = gsl_matrix_int_alloc(1, nx);
ind     = mxGetPr(prhs[0]);
for(j = 0; j < nx; j++){
gsl_matrix_int_set(index, 0, (const size_t)j, (int)ind[j]);
}

// Retrieve matrices Ai and Bi

for(k = 0; k < nrhs-4; k++){
  

if(! mxIsEmpty(prhs[k+1]))
{
  /* Get the length of each input vector. */
  mx = mxGetM(prhs[k+1]);
  nx = mxGetN(prhs[k+1]);
  /* Check input parameters size */

if(k == 0 && nx != dim){mexErrMsgTxt("Missing Fox H argument(s) z\n");}
if(k > 0 && nx ==1){
         mexPrintf("Error in parameter # %d\n",k+2); 
         mexErrMsgTxt("Input size is incorrect\n");            
         return;
         }

if(mxIsComplex(prhs[k+1])){
  /* Get pointers to real and imaginary parts of the inputs. */
  xr = mxGetPr(prhs[k+1]);
  xi = mxGetPi(prhs[k+1]);

  Arg[k] = gsl_matrix_complex_alloc(mx, nx);

for(i = 0; i < mx; i++)
  {
      for(j = 0; j < nx; j++)
      {
          GSL_SET_COMPLEX(x,xr[i + mx*j],xi[i + mx*j]);
          gsl_matrix_complex_set(Arg[k], 
                     (const size_t)i, (const size_t)j, x);
      }
  }

}//end of if

else{

/* Get pointers to real part */
  xr = mxGetPr(prhs[k+1]);

  Arg[k] = gsl_matrix_complex_alloc(mx, nx);

for(i = 0; i < mx; i++)
  {
      for(j = 0; j < nx; j++)
      {
          GSL_SET_COMPLEX(x,xr[i + mx*j],0.0);
          gsl_matrix_complex_set(Arg[k], 
                     (const size_t)i, (const size_t)j, x);
      }
  }

}//end of else
}

else{
  Arg[k] = gsl_matrix_complex_alloc(2,1);
  gsl_matrix_complex_set_all(Arg[k],GSL_COMPLEX_ONE);
}

}// end of for(k=0...

// Retrieve integration intervals
  if(mxGetN(prhs[nrhs-3]) < dim || mxGetM(prhs[nrhs-3]) < 2){
                  mexErrMsgTxt("Contour matrix size incorrect"); 
                  }
  xl = gsl_vector_complex_alloc(dim);
  xu = gsl_vector_complex_alloc(dim);

/* Get pointers to real and imaginary parts of the inputs. */
  xr = mxGetPr(prhs[nrhs-3]);
  xi = mxGetPi(prhs[nrhs-3]);

//Initialize integration domains
for (i = 0; i < dim; i++) {
  gsl_vector_complex_set (xl, i, gsl_complex_rect(xr[2*i],xi[2*i]));
  gsl_vector_complex_set (xu, i, gsl_complex_rect(xr[2*i+1],xi[2*i+1]));
}

// Retrive max_call and tolerence
if(!mxIsComplex(prhs[nrhs-2]) && !mxIsComplex(prhs[nrhs-1])){
   max_call = mxGetPr(prhs[nrhs-2]);
   tol      = mxGetPr(prhs[nrhs-1]);
  }
else{mexErrMsgTxt("MaxFunEval and AbsTol must be real");}
  
//Perform integration
  gsl_complex result = GSL_COMPLEX_ZERO, error= GSL_COMPLEX_ZERO;
  gsl_qrng* qrng = gsl_qrng_alloc(gsl_qrng_reversehalton, dim);
/*Available sequences: gsl_qrng_halton, gsl_qrng_sobol, gsl_qrng_niederreiter_2, gsl_qrng_reversehalton*/
  quasi_monte_state* s = quasi_monte_alloc(dim);
  quasi_monte_integrate(xl, xu, dim, max_call[0], 0, tol[0], qrng, s, &result, &error, index, Arg);
  quasi_monte_free(s);
  gsl_qrng_free(qrng);

/* Create a new complex array and set the output pointer to it. */
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
  zr = mxGetPr(plhs[0]);
  zi = mxGetPi(plhs[0]);
  zr[0] = GSL_REAL(result);
  zi[0] = GSL_IMAG(result);
}


