#ifndef ADDITIONAL_MATH_
#define ADDITIONAL_MATH_
#include <math.h>
#include "definitions.h"


/* Function prototypes for linear assembly routines*/
extern float mag_sqrd(COMPLEX z);
extern COMPLEX cmplx_mult_add(COMPLEX w, COMPLEX x, COMPLEX y, COMPLEX z, float *real, float *imag);

/* Function prototypes additional_math.c functions*/
float FastInv4over3(float x);
float FastInvSqrt(float x);
void eig(COMPLEX *M, COMPLEX *eigvals, COMPLEX *eigvecs);
float mag(COMPLEX x);
double mag_dbl(COMPLEX_DBL x);
float arg(COMPLEX x);
COMPLEX conj(COMPLEX z);
COMPLEX_DBL conj_dbl(COMPLEX_DBL z);
void inv_2x2(COMPLEX *A, COMPLEX *out);
COMPLEX cmplx_mult(COMPLEX z1, COMPLEX z2);
COMPLEX_DBL cmplx_mult_dbl(COMPLEX_DBL z1, COMPLEX_DBL z2);
COMPLEX cmplx_add(COMPLEX z1, COMPLEX z2);
COMPLEX cmplx_minus(COMPLEX z1, COMPLEX z2);
void eig_dbl(COMPLEX_DBL *M, COMPLEX_DBL *eigvals, COMPLEX_DBL *eigvecs);


#endif /*ADDITIONAL_MATH_*/
