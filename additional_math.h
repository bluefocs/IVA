#ifndef ADDITIONAL_MATH_
#define ADDITIONAL_MATH_
#include <math.h>
#include "definitions.h"


/* Function prototypes*/
void eig(COMPLEX *M, COMPLEX *eigvals, COMPLEX *eigvecs);
float mag(COMPLEX x);
float arg(COMPLEX x);
COMPLEX conj(COMPLEX z);
void inv_2x2(COMPLEX *A, COMPLEX *out);
COMPLEX cmplx_mult(COMPLEX z1, COMPLEX z2);
COMPLEX cmplx_add(COMPLEX z1, COMPLEX z2);
COMPLEX cmplx_minus(COMPLEX z1, COMPLEX z2);


#endif /*ADDITIONAL_MATH_*/
