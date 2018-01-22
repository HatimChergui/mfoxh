/*
 * Copyright 2018 Hatim Chergui, Mustapha Benjillali, 
 * and Mohamed-Slim Alouini. Contact email <chergui@ieee.org>
 *
 * If you use this software or any (modified) part of it, please cite 
 * the submitted paper: "K-Factor-Based Analysis of XLOS Service 
 * Probability in 5G Ultra-Dense Networks". [Online]. Available: arxiv.org/
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

#ifndef _MFOX_H_INCLUDE
#define _MFOX_H_INCLUDE

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_matrix_int.h>
#include <gsl/gsl_qrng.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_nan.h>
//#include "gsl_complex_monte.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  size_t dim;
  double *x;
} quasi_monte_state;

quasi_monte_state* quasi_monte_alloc(size_t dim);
int quasi_monte_init(quasi_monte_state* s);
int quasi_monte_integrate(gsl_vector_complex *xl, gsl_vector_complex *xu, size_t dim, size_t max_calls, double max_relerr, double max_abserr, gsl_qrng* r, quasi_monte_state* state, gsl_complex *result, gsl_complex *abserr, gsl_matrix_int *index, gsl_matrix_complex *Arg[20]);
void quasi_monte_free(quasi_monte_state* s);
int gammaz(gsl_complex z, gsl_complex* b);
int ProdArg(gsl_matrix_complex *A, const gsl_vector_complex *x, gsl_complex *z);
int Phi(gsl_matrix_complex *A, gsl_matrix_complex *B, const gsl_complex x, size_t m, size_t n, gsl_complex *z);
int Psi(gsl_matrix_complex *A, gsl_matrix_complex *B, const gsl_vector_complex *x, size_t n, gsl_complex *z);
gsl_complex mFoxIntegrand(const gsl_vector_complex *x, size_t dim, gsl_matrix_int *index, gsl_matrix_complex *Arg[20]);

#ifdef __cplusplus
}
#endif

#endif

