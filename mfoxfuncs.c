/*
 * Copyright 2018 Hatim Chergui, Mustapha Benjillali, 
 * and Mohamed-Slim Alouini. Contact email <chergui@ieee.org>
 ****************************************************************************************************************
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
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_matrix_int.h>
#include "mfox.h"

//###################### Gammaz Function ####################
#define GSL_FOXH_FACTOR (gsl_complex_inverse(gsl_complex_rect(0.0,2*M_PI)))

int gammaz(gsl_complex z, gsl_complex* b)
{
  gsl_sf_result lnr, arg;
  gsl_sf_lngamma_complex_e (GSL_REAL(z), GSL_IMAG(z), &lnr, &arg);
  *b = gsl_complex_polar(exp(lnr.val), arg.val);
  return GSL_SUCCESS;
  
}

//##################### ProdArg Function ####################
int ProdArg(gsl_matrix_complex *A, const gsl_vector_complex *x, gsl_complex *z){

size_t j, dim = A->size2;
gsl_complex zt;
*z = GSL_COMPLEX_ONE;

for(j = 0; j< dim; j++){

zt = gsl_matrix_complex_get(A, 0, j);

*z = gsl_complex_mul(*z, gsl_complex_pow(zt, gsl_vector_complex_get(x,j)));
}
}

//###################### Phi Function ####################
int Phi(gsl_matrix_complex *A, gsl_matrix_complex *B, const gsl_complex x, size_t m, size_t n, gsl_complex *z)
{

size_t p, q, dim, i, j;
gsl_complex zn = GSL_COMPLEX_ONE, zp = GSL_COMPLEX_ONE, zm = GSL_COMPLEX_ONE, zq = GSL_COMPLEX_ONE, tmpn = GSL_COMPLEX_ONE, tmpp = GSL_COMPLEX_ONE, tmpm = GSL_COMPLEX_ONE, tmpq = GSL_COMPLEX_ONE;
gsl_vector_complex *In = NULL, *Iq = NULL;
gsl_vector_complex_view an, An, ap, Ap, bm, Bm, bq, Bq;
gsl_matrix_complex * E, *F;
E = gsl_matrix_complex_alloc(A->size1,A->size2);
F = gsl_matrix_complex_alloc(B->size1,B->size2);
gsl_matrix_complex_memcpy(E, A);
gsl_matrix_complex_memcpy(F, B);

p   = E->size1;
q   = F->size1;

if(E->size2 == 1){zn = GSL_COMPLEX_ONE; zp = GSL_COMPLEX_ONE; }
if(F->size2 == 1){zm = GSL_COMPLEX_ONE; zq = GSL_COMPLEX_ONE; }

if(E->size2 != 1){

if(n > p){mexErrMsgTxt("Index ni should be <= pi");}

if (n == 0){zn = GSL_COMPLEX_ONE;}

if(n > 0){ 
           In = gsl_vector_complex_alloc(n);
           gsl_vector_complex_set_all(In, GSL_COMPLEX_ONE);
           an = gsl_matrix_complex_subcolumn (E, 0, 0, n);
           An = gsl_matrix_complex_subcolumn (E, 1, 0, n);
           gsl_vector_complex_scale(&An.vector, x);
           gsl_vector_complex_sub(In, &an.vector);
           gsl_vector_complex_add(In, &An.vector);

           for(i = 0; i < n; i++){
            gammaz(gsl_vector_complex_get(In,i), &tmpn);
            zn = gsl_complex_mul(zn, tmpn) ;
           }
}

if (p == n || p ==0){zp = GSL_COMPLEX_ONE;}

if(n < p){
           ap = gsl_matrix_complex_subcolumn (E, 0, n, p-n);
           Ap = gsl_matrix_complex_subcolumn (E, 1, n, p-n);
           gsl_vector_complex_scale(&Ap.vector, x);
           gsl_vector_complex_sub(&ap.vector, &Ap.vector);

           for(i = 0; i < p-n; i++){
             gammaz(gsl_vector_complex_get(&ap.vector,i), &tmpp);
             zp = gsl_complex_mul(zp, tmpp) ;
           }
}
}

if(F->size2 != 1){

if(m > q){mexErrMsgTxt("Index mi should be <= qi");}

if (m == 0){zm = GSL_COMPLEX_ONE;}

if(m > 0){
           bm = gsl_matrix_complex_subcolumn (F, 0, 0, m);
           Bm = gsl_matrix_complex_subcolumn (F, 1, 0, m);
           gsl_vector_complex_scale(&Bm.vector, x);
           gsl_vector_complex_sub(&bm.vector,&Bm.vector);

           for(i=0; i < m; i++){
            gammaz(gsl_vector_complex_get(&bm.vector,i), &tmpm);
            zm = gsl_complex_mul(zm, tmpm);
           }
}

if (q == m || q ==0){zq = GSL_COMPLEX_ONE;}

if(m < q){
           Iq = gsl_vector_complex_alloc(q-m);
           gsl_vector_complex_set_all(Iq, GSL_COMPLEX_ONE);
           bq = gsl_matrix_complex_subcolumn (F, 0, m, q-m);
           Bq = gsl_matrix_complex_subcolumn (F, 1, m, q-m);
           gsl_vector_complex_scale(&Bq.vector, x);
           gsl_vector_complex_sub(Iq, &bq.vector);
           gsl_vector_complex_add(Iq, &Bq.vector);

           for(i=0; i < q-m; i++){
            gammaz(gsl_vector_complex_get(Iq,i), &tmpq);
            zq = gsl_complex_mul(zq, tmpq) ;
           }
}
}

*z = gsl_complex_div(gsl_complex_mul(zn,zm), gsl_complex_mul(zp,zq));
return(EXIT_SUCCESS);
}



//##################### Psi Function ######################

int Psi(gsl_matrix_complex *A, gsl_matrix_complex *B, const gsl_vector_complex *x, size_t n, gsl_complex *z)
{

size_t p, q, dim, i, j;
gsl_complex zn = GSL_COMPLEX_ONE, zp = GSL_COMPLEX_ONE, zm = GSL_COMPLEX_ONE, zq = GSL_COMPLEX_ONE, tn = GSL_COMPLEX_ONE, tp = GSL_COMPLEX_ONE, tq = GSL_COMPLEX_ONE;
gsl_vector_complex *In = NULL, *Iq = NULL;
gsl_matrix_complex *E, *F;
E = gsl_matrix_complex_alloc(A->size1,A->size2);
F = gsl_matrix_complex_alloc(B->size1,B->size2);
gsl_matrix_complex_memcpy(E, A);
gsl_matrix_complex_memcpy(F, B);
gsl_vector_complex_view an, ap, bq;
gsl_matrix_complex_view An, Ap, Bq;

p   = E->size1;
q   = F->size1;
dim = E->size2 - 1;

if(E->size2 == 1){zn = GSL_COMPLEX_ONE; zp = GSL_COMPLEX_ONE; }

if(F->size2 == 1){zq = GSL_COMPLEX_ONE; }

if(E->size2 != 1){
if(n > p){mexErrMsgTxt("Index ni should be <= pi");}
if(n == 0){zn = GSL_COMPLEX_ONE;}
if(n > 0){
          In = gsl_vector_complex_alloc(n);
          gsl_vector_complex_set_all(In, GSL_COMPLEX_ONE);
          an = gsl_matrix_complex_subcolumn (E, 0, 0, n);
          An = gsl_matrix_complex_submatrix (E, 0, 1, n, dim);
          gsl_vector_complex_sub(In, &an.vector);
          gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, &An.matrix, x, GSL_COMPLEX_ONE, In);
          for(i = 0; i < n; i++){
           gammaz(gsl_vector_complex_get(In,i), &tn);
           zn = gsl_complex_mul(zn, tn);}
}

if (p == n || p ==0){zp = GSL_COMPLEX_ONE;}      
if(n < p){
          ap = gsl_matrix_complex_subcolumn (E, 0, n, p-n);
          Ap = gsl_matrix_complex_submatrix (E, n, 1, p-n, dim);
          gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_NEGONE, &Ap.matrix, x, GSL_COMPLEX_ONE, &ap.vector);
          for(i = 0; i < p-n; i++){
           gammaz(gsl_vector_complex_get(&ap.vector,i), &tp);
           zp = gsl_complex_mul(zp, tp);
          }

}

}

if(F->size2 != 1){
          Iq = gsl_vector_complex_alloc(q);
          gsl_vector_complex_set_all(Iq, GSL_COMPLEX_ONE);
          bq = gsl_matrix_complex_subcolumn (F, 0, 0, q);
          Bq = gsl_matrix_complex_submatrix (F, 0, 1, q, dim);
          gsl_vector_complex_sub(Iq, &bq.vector);
          gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, &Bq.matrix, x, GSL_COMPLEX_ONE, Iq);
          for(i=0; i < q; i++){
            gammaz(gsl_vector_complex_get(Iq,i), &tq);
            zq = gsl_complex_mul(zq, tq) ;
          }
}

*z = gsl_complex_div(zn, gsl_complex_mul(zp,zq));
return(EXIT_SUCCESS);

}


//###################### mFoxIntegrand #####################
gsl_complex mFoxIntegrand(const gsl_vector_complex *x, size_t dim, gsl_matrix_int *index, gsl_matrix_complex *Arg[20]) {

size_t i, mx, nx, narg = 2*dim+3;
gsl_complex z, zarg, zi, zphi = GSL_COMPLEX_ONE;
   
// Parameters [0, n, m1, n1, ..., mM, nM] are stored in matrix index

// Call ProdArg function
ProdArg(Arg[0], x, &zarg);

// Call Psi function
Psi(Arg[1], Arg[2], x, (size_t) gsl_matrix_int_get(index,0,1), &z);

// Call Phi function
for(i = 0; i < dim; i++){
Phi(Arg[2*i+3], Arg[2*i+4], gsl_vector_complex_get(x,i), (size_t) gsl_matrix_int_get(index,0,2*i+2), (size_t) gsl_matrix_int_get(index,0,2*i+3), &zi);
zphi = gsl_complex_mul(zphi,zi);
}

gsl_complex NFOXH_Factor = gsl_complex_pow(GSL_FOXH_FACTOR, gsl_complex_rect((double)dim,0.0));

return gsl_complex_mul(gsl_complex_mul(gsl_complex_mul(z,zphi), zarg), NFOXH_Factor);

}

//####################### Complex QMC ########################

quasi_monte_state* quasi_monte_alloc(size_t dim) {
    quasi_monte_state* s = (quasi_monte_state*)malloc(sizeof(quasi_monte_state));
    if (s == NULL) {
        GSL_ERROR_VAL("failed to allocate space for quasi_monte_state", GSL_ENOMEM, 0);
    }
    s->x = (double*)malloc(dim * sizeof(double));
    if (s->x == NULL) {
        free(s);
        GSL_ERROR_VAL("failed to allocate space for working vector", GSL_ENOMEM, 0);
    }
    s->dim = dim;
    return s;
}

int quasi_monte_init(quasi_monte_state* s) {
    size_t i;
    for (i = 0; i < s->dim; i++) {
        s->x[i] = 0;
    }
    return GSL_SUCCESS;
}

int quasi_monte_integrate(gsl_vector_complex *xl, gsl_vector_complex *xu, size_t dim, size_t max_calls, double max_relerr, double max_abserr, gsl_qrng* r, quasi_monte_state* state, gsl_complex *result, gsl_complex *abserr, gsl_matrix_int *index, gsl_matrix_complex *Arg[20]){

gsl_complex volume = GSL_COMPLEX_ONE, mean = GSL_COMPLEX_ZERO, variance_sum = GSL_COMPLEX_ZERO, running_avg_abserr = gsl_complex_rect(GSL_POSINF,GSL_POSINF);
//double *xlr, *xli, *xur, *xui; 
size_t n, i;
gsl_vector_complex *y = NULL;
y = gsl_vector_complex_alloc(dim);
gsl_vector_complex_set_all(y, GSL_COMPLEX_ZERO);
double* x = state->x;
gsl_complex fval;
    

    // Check that the dimensionalities match
    if (dim != state->dim) {
        char errmsg[80];
        snprintf(errmsg, 80, "number of dimensions %ud doesn't match allocated size %ud", (unsigned int)dim, (unsigned int)(state->dim));
        GSL_ERROR(errmsg, GSL_EINVAL);
    }

    // Check the bounds for validity

    for (i = 0; i < dim; i++) {
        gsl_complex tl = gsl_vector_complex_get(xl,i);
        gsl_complex tu = gsl_vector_complex_get(xu,i);

if (GSL_REAL(tu) < GSL_REAL(tl) || GSL_IMAG(tu) <= GSL_IMAG(tl)) {                           
            char errmsg[80];
            snprintf(errmsg, 80, "lower limit not less than upper limit");
            GSL_ERROR(errmsg, GSL_EINVAL);
        }
        if (GSL_REAL(tu) - GSL_REAL(tl) > GSL_DBL_MAX || GSL_IMAG(tu) - GSL_IMAG(tl) > GSL_DBL_MAX) {
            char errmsg[80];
            snprintf(errmsg, 80, "integration range is larger than limit %f", GSL_DBL_MAX);
            GSL_ERROR(errmsg, GSL_EINVAL);
        }
    }
    // Compute the volume of the region
    for (i = 0; i < dim; i++) {
gsl_complex vu = gsl_vector_complex_get(xu,i);
gsl_complex vl = gsl_vector_complex_get(xl,i);

        volume = gsl_complex_mul(volume, gsl_complex_sub(vu,vl));
    }
    for (n = 0; n < max_calls; n++) {
        // Choose a quasirandom point in the integration region
        gsl_qrng_get(r, x);
        for (i = 0; i < dim; i++) {
gsl_complex vu = gsl_vector_complex_get(xu,i);
gsl_complex vl = gsl_vector_complex_get(xl,i);
gsl_vector_complex_set(y, i, gsl_complex_add(vl , gsl_complex_mul_real (gsl_complex_sub(vu,vl), x[i])));
        }

        {
           fval = mFoxIntegrand(y, dim, index, Arg);
           gsl_complex d = gsl_complex_sub(fval , mean);
           mean = gsl_complex_add(mean, gsl_complex_div_real(d ,(double)(n + 1.0)));
           variance_sum = gsl_complex_add(variance_sum,gsl_complex_mul_real(gsl_complex_mul(d,d), (double)(n/(n + 1.0))));
        }
        if (n > 1) {
            running_avg_abserr = gsl_complex_sqrt(gsl_complex_div_real(variance_sum,(double)((n + 1.0) * n)));

 gsl_complex tol = gsl_complex_mul(volume, running_avg_abserr );


            if (fabs(GSL_REAL(tol)) < max_abserr && fabs(GSL_IMAG(tol)) < max_abserr) {
                break;
            }
            /*if (gsl_complex_abs(gsl_complex_div(running_avg_abserr, mean)) < max_relerr) {
                break;
            }*/
        }
    }
    *result = gsl_complex_mul(volume , mean);
    *abserr = gsl_complex_mul(volume , running_avg_abserr);
    return GSL_SUCCESS;
}

void quasi_monte_free(quasi_monte_state* s) {
    if (s == NULL) {
        return;
    }
    free(s->x);
    free(s);
}


