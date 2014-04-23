#include "additional_math.h"

float FastInv4over3(float x) 
{
	//float xhalf = 0.5f * x;
	int i = *(int*)&x;
	i = 0x941D3679 - (1.333f * i);
	x = *(float*)&i;
	//x = x*(1.5f-(xhalf*x*x)); // Can't do Newton-Raphson with 4/3 as defeats point of the algorithm
	return x;
}

float FastInvSqrt(float x) 
{
	float xhalf = 0.5f * x;
	int i = *(int*)&x;         // evil floating point bit level hacking
	i = 0x5f3759df - (i >> 1);  // what the fuck?
	x = *(float*)&i;
	x = x*(1.5f-(xhalf*x*x));
	return x;
}

void inv_2x2(COMPLEX *A, COMPLEX *out)
{
	// Inverse of complex valued 2x2 matrix - assume length of A is 4
	COMPLEX det,detinv; // Determinate
	
	det = cmplx_minus(cmplx_mult(A[0],A[3]),cmplx_mult(A[1],A[2]));
	detinv.real = (1/mag(det)) * cos(arg(det));
	detinv.imag = (1/mag(det)) * sin(arg(det)) * (-1);
	
	out[0].real = A[3].real * detinv.real;
	out[0].imag = A[3].imag * detinv.imag;
	out[3].real = A[0].real * detinv.real;
	out[3].imag = A[0].imag * detinv.imag;
	out[1].real = A[1].real * detinv.real * (-1.0);
	out[1].imag = A[1].imag * detinv.imag * (-1.0);
	out[2].real = A[2].real * detinv.real * (-1.0);
	out[2].imag = A[2].imag * detinv.imag * (-1.0);
}


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
	//r = pow((double)r, 0.5);
	r = sqrt(r);
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
	
	/*
	//Normalise
	norm_fact = sqrt(mag(eigvecs[0])*mag(eigvecs[0]) + mag(eigvecs[2])*mag(eigvecs[2]));
	eigvecs[0].real = eigvecs[0].real / norm_fact;
	eigvecs[0].imag = eigvecs[0].imag / norm_fact;
	eigvecs[2].real = eigvecs[2].real / norm_fact;
	eigvecs[2].imag = eigvecs[2].imag / norm_fact;
	
	norm_fact = sqrt(mag(eigvecs[1])*mag(eigvecs[1]) + mag(eigvecs[3])*mag(eigvecs[3]));
	eigvecs[1].real = eigvecs[1].real / norm_fact;
	eigvecs[1].imag = eigvecs[1].imag / norm_fact;
	eigvecs[3].real = eigvecs[3].real / norm_fact;
	eigvecs[3].imag = eigvecs[3].imag / norm_fact;
	*/
	
	//Normalise - may not be as accurate as the above method
	norm_fact = mag(eigvecs[0])*mag(eigvecs[0]) + mag(eigvecs[2])*mag(eigvecs[2]);
	eigvecs[0].real = eigvecs[0].real * FastInvSqrt( norm_fact );
	eigvecs[0].imag = eigvecs[0].imag * FastInvSqrt( norm_fact );
	eigvecs[2].real = eigvecs[2].real * FastInvSqrt( norm_fact );
	eigvecs[2].imag = eigvecs[2].imag * FastInvSqrt( norm_fact );
	
	norm_fact = mag(eigvecs[1])*mag(eigvecs[1]) + mag(eigvecs[3])*mag(eigvecs[3]);
	eigvecs[1].real = eigvecs[1].real * FastInvSqrt( norm_fact );
	eigvecs[1].imag = eigvecs[1].imag * FastInvSqrt( norm_fact );
	eigvecs[3].real = eigvecs[3].real * FastInvSqrt( norm_fact );
	eigvecs[3].imag = eigvecs[3].imag * FastInvSqrt( norm_fact );
}
float mag(COMPLEX x)
{
	return sqrt(x.real*x.real + x.imag*x.imag);
}
float arg(COMPLEX x)
{
	return atan(x.imag/x.real);
}
COMPLEX conj(COMPLEX z)
{
	COMPLEX out;
	out.real = z.real;
	out.imag = (-1)*z.imag;
	return out;
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
/*
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
*/
