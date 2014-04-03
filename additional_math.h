#ifndef ADDITIONAL_MATH_
#define ADDITIONAL_MATH_
#include <math.h>
/* Function prototypes*/
void eig(COMPLEX *X, COMPLEX *eigvals, COMPLEX *eigvecs);
float mag(COMPLEX x);
float arg(COMPLEX x);
COMPLEX cmplx_mult(COMPLEX z1, COMPLEX z2);
COMPLEX cmplx_add(COMPLEX z1, COMPLEX z2);
COMPLEX cmplx_minus(COMPLEX z1, COMPLEX z2);

void eig(COMPLEX *M, COMPLEX *eigvals, COMPLEX *eigvecs)
{
	// Calculates the eigenvalues and eigenvectors of 2*2 matrix
	// Using the equations here: http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
	// Accepts complex and real numbers, ie. complex part is set to zero.
	
	//COMPLEX D;
	//COMPLEX T;
	
	double Dr = 0.0, Di = 0.0; // Real and imaginary parts as doubles
	double Tr = 0.0, Ti = 0.0; // .. necessary to have double precision here.
	
	float norm_fact, arg;
	double r=0.0; // Double so that it doesn't go to inf
	float realpart=0.0, imagpart=0.0;
	
	/* Work out determinate (2x2) */
	Dr = ((double)(M[0].real)*(double)(M[3].real)) - ((double)(M[0].imag)*(double)(M[3].imag)) - ((double)(M[1].real)*(double)(M[2].real)) + ((double)(M[1].imag)*(double)(M[2].imag));  
	Di = ((double)(M[0].real)*(double)(M[3].imag)) + ((double)(M[0].imag)*(double)(M[3].real)) - ((double)(M[2].real)*(double)(M[1].imag)) - ((double)(M[2].imag)*(double)(M[1].real));  
	
	Tr = M[0].real + M[3].real; //Trace real part
	Ti = M[0].imag + M[3].imag;
	
// for the eigenvalues to be found need to express the bit under the square
// root as (stuff) + j(stuff) so that it can be rewritten sqrt(r)(cos(0.5*arg) + sin(0.5*arg))
	
//	realpart = (float)((T.real*T.real)-(T.imag*T.imag))/4.0;
//	realpart -= D.real; 
//	imagpart = (0.5*T.imag*T.real) - D.imag;
	
	realpart = (float)((Tr*Tr)-(Ti*Ti))/4.0;
	realpart -= Dr; 
	imagpart = (0.5*Ti*Tr) - Di;
	
	
	arg = atan(imagpart/realpart);
	r = ((double)realpart*(double)realpart) + ((double)imagpart*(double)imagpart);
	//arg = tan((0.5*T.imag*T.real - D.imag)/((pow(T.real,2)-pow(T.imag,2))/4 - D.real));
	//r = pow((pow(T.real,2)-pow(T.imag,2))/4 - D.real,2) + pow((0.5*T.imag*T.real - D.imag),2);
	//r = sqrt((double)r);
	r = pow((double)r, 0.5);
	eigvals[0].real = (Tr/2.0) + sqrt(r)*cos(0.5*arg);
	eigvals[0].imag = (Ti/2.0) + sqrt(r)*sin(0.5*arg);
	eigvals[1].real = (Tr/2.0) - sqrt(r)*cos(0.5*arg);
	eigvals[1].imag = (Ti/2.0) - sqrt(r)*sin(0.5*arg);
///	eigvals[0] = (T/2) + sqrt((pow(T,2)/4) - D);//Real case
//	eigvals[1] = (T/2) - sqrt((pow(T,2)/4) - D);//Real case
	
	
	if((M[2].real==0.0) && (M[1].real==0.0) && (M[2].imag==0.0) && (M[1].imag==0.0))
	{
		eigvecs[0].real=1; // 2*2 identity matrix
		eigvecs[0].imag=0;
		eigvecs[1].real=0;
		eigvecs[1].imag=0;
		eigvecs[2].real=0;
		eigvecs[2].imag=0;
		eigvecs[3].real=1;
		eigvecs[3].imag=0;
	}
	else if((M[2].real!=0.0) && (M[2].imag!=0.0))// if c is not zero
	{
		eigvecs[0].real = eigvals[0].real-M[3].real;
		eigvecs[0].imag = eigvals[0].imag-M[3].imag;
		eigvecs[1].real = eigvals[1].real-M[3].real;
		eigvecs[1].imag = eigvals[1].imag-M[3].imag;
		eigvecs[2].real = M[2].real;//c
		eigvecs[2].imag = M[2].imag;//c
		eigvecs[3].real = M[2].real;//c
		eigvecs[3].imag = M[2].imag;//c
	}
	else // if b!=0
	{
		eigvecs[0].real = M[1].real;//b
		eigvecs[0].imag = M[1].imag;//b
		eigvecs[1].real = M[1].real;//b
		eigvecs[1].imag = M[1].imag;//b
		eigvecs[2].real = eigvals[0].real-M[0].real;
		eigvecs[2].imag = eigvals[0].imag-M[0].imag;
		eigvecs[3].real = eigvals[1].real-M[0].real;
		eigvecs[3].imag = eigvals[1].imag-M[0].imag;
	}
	
	
	//Normalise
//	norm_fact = sqrt(pow(mag(eigvecs[0]),(double)2.0) + pow(mag(eigvecs[2]),(double)2.0));// normalisation factor	
	norm_fact = sqrt(mag(eigvecs[0])*mag(eigvecs[0]) + mag(eigvecs[2])*mag(eigvecs[2]));
	eigvecs[0].real = eigvecs[0].real / norm_fact;
	eigvecs[0].imag = eigvecs[0].imag / norm_fact;
	eigvecs[2].real = eigvecs[2].real / norm_fact;
	eigvecs[2].imag = eigvecs[2].imag / norm_fact;
	
//	norm_fact = sqrt(pow(mag(eigvecs[1]),(double)2.0) + pow(mag(eigvecs[3]),(double)2.0));// normalisation factor	
	norm_fact = sqrt(mag(eigvecs[1])*mag(eigvecs[1]) + mag(eigvecs[3])*mag(eigvecs[3]));
	eigvecs[1].real = eigvecs[1].real / norm_fact;
	eigvecs[1].imag = eigvecs[1].imag / norm_fact;
	eigvecs[3].real = eigvecs[3].real / norm_fact;
	eigvecs[3].imag = eigvecs[3].imag / norm_fact;
}
float mag(COMPLEX x)
{
	return sqrt(x.real*x.real + x.imag*x.imag);
}
float arg(COMPLEX x)
{
	return atan(x.imag/x.real);
}
COMPLEX cmplx_mult(COMPLEX z1, COMPLEX z2)
{
	COMPLEX out;
	float aminusb;
	//aminusb = z1.real + z2.real;// Why on earth was I doing this????
	aminusb = (z1.real - z1.imag) * z2.imag;
	// Uses an identity to reduce the number of multiplications
	out.real = z1.real*(z2.real - z2.imag) + aminusb;
	out.imag = z1.imag*(z2.real + z2.imag) + aminusb;
	//out.real = z1.real*z2.real - (z1.imag*z2.imag);
	//out.imag = z1.real*z2.imag + (z1.imag*z2.real);
	
	return out;
}
COMPLEX cmplx_add(COMPLEX z1, COMPLEX z2)
{
	COMPLEX out;
	
	out.real = z1.real+z2.real;// Real
	out.imag = z1.imag+z2.imag; // Imaginary
	
	return out;
}
COMPLEX cmplx_minus(COMPLEX z1, COMPLEX z2)
{
	COMPLEX out;
	
	out.real = z1.real-z2.real;// Real
	out.imag = z1.imag-z2.imag; // Imaginary
	
	return out;
}
void COMPLEX_sp_mat_mul(COMPLEX *x, int r1, int c1, COMPLEX *y, int c2, COMPLEX *r)
{
	// Similar to DSPF_sp_mat_mul, but for the complex case
	short i, j, k;
	COMPLEX sum;
	// Multiply each row in x by each column in y.
	// The product of row m in x and column n in y is placed
	// in position (m,n) in the result.
	for (i = 0; i < r1; i++)
	for (j = 0; j < c2; j++)
	{
		sum.real = 0.0;
		sum.imag = 0.0;
		for (k = 0; k < c1; k++)
		{		
//			sum += x[k + i*c1] * y[j + k*c2];
			sum = cmplx_add(sum,cmplx_mult(x[k + i*c1], y[j + k*c2]));
		}
		r[j + i*c2] = sum;
	}
}
void COMPLEX_sp_mat_trans(COMPLEX *x, int rows, int cols, COMPLEX *r, unsigned char herm_flag)
{
	short i,j;
	for(i=0; i<cols; i++)
	for(j=0; j<rows; j++)
	{
		if(herm_flag==1)// Heritian martix flag is set
		{
			r[i * rows + j].real = x[i + cols * j].real;
			r[i * rows + j].imag = x[i + cols * j].imag * (-1);
		}
		else
		{
			r[i * rows + j].real = x[i + cols * j].real;
			r[i * rows + j].imag = x[i + cols * j].imag;
		}
	}
}


#endif /*ADDITIONAL_MATH_*/