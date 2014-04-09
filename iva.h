#ifndef IVA_H_
#define IVA_H_
#include "definitions.h"

COMPLEX S[N2 * 2 * TIME_BLOCKS_50PC]; //

void iva(COMPLEX *Xp, COMPLEX *Wp, unsigned short nfreq)
{
	short k=0, m=0;
	unsigned short maxiter=1000, iter=0;
	float mu=0.1;
	COMPLEX detWp;
	float Ssq[TIME_BLOCKS_50PC * NSOURCES];
	double epsilon = 0.000001;
	COMPLEX W_temp[4], W_new[4];//, detWp;	
	COMPLEX Phi[TIME_BLOCKS_50PC * NSOURCES];
	double SumSsq=0.0, dObj=0.0, pObj=0.0, Obj=0.0, dlw=0.0, tol = 0.000001, comparison=0.0;
	
	for(iter=0;iter<maxiter;iter++)
	{
		dlw = 0;
		
		for(k=0;k<nfreq;k++)
		{
			//COMPLEX_sp_mat_mul(&(Wp[4*k + 0]), nsou, nsou, &([CH1 + N2*k + m]), nblocks, &(S[k][0][0]));	
			for(m=0;m<TIME_BLOCKS_50PC;m++)// 2 by many matrix multiplied by many by 2 matrix
			{
				S[CH1 + N2*m + k] = cmplx_add(cmplx_mult(Wp[4*k + 0], Xp[CH1 + N2*m + k]), cmplx_mult(Wp[4*k + 1], Xp[CH2 + N2*m + k]));
				S[CH2 + N2*m + k] = cmplx_add(cmplx_mult(Wp[4*k + 2], Xp[CH1 + N2*m + k]), cmplx_mult(Wp[4*k + 3], Xp[CH2 + N2*m + k]));
			}
		}
		
		SumSsq=0.0;
		// Calculate non-linear function?
		for(m=0;m<TIME_BLOCKS_50PC;m++) // Summnation loop - Can this be sped up ? and be done with double precision?
		{
			for(k=0;k<nfreq;k++)
			{
				Ssq[ m ] += pow(mag(S[CH1 + N2*m + k]), 2.0);
				Ssq[TIME_BLOCKS_50PC+m] += pow(mag(S[CH2 + N2*m + k]), 2.0);
			}
			
			Ssq[ m ] = sqrt(Ssq[ m ]); // In the future change ^0.5 to ^0.666. Important line! 
			Ssq[TIME_BLOCKS_50PC+m] = sqrt(Ssq[TIME_BLOCKS_50PC+m]); // Channel 2
			
			// Calculate the sum of all the values in Ssq before the inverse is taken, this is used in the break condition below
			SumSsq += Ssq[ m ] + Ssq[TIME_BLOCKS_50PC+m];
		}
		for(m=0;m<TIME_BLOCKS_50PC;m++) // Take the summnation of 
		{
			Ssq[ m ] = 1.0/(Ssq[ m ] + epsilon);	// Channel 1 - Ssq1 in MATLAB code
			Ssq[TIME_BLOCKS_50PC+m]=1.0/(Ssq[TIME_BLOCKS_50PC+m] + epsilon);	//Channel 2
		}
		
		

		for(k=0;k<nfreq;k++)
		{
			//Calculate multivariate score function and gradients
			for(m=0; m<TIME_BLOCKS_50PC; m++)
			{
				// Unrolled the source loop
					
				Phi[ CH1 + m ].real = S[CH1 + N2*m + k].real * Ssq[CH1 + m]; // Phi exists at each frequency bin for each channel
				Phi[ CH1 + m ].imag = S[CH1 + N2*m + k].imag * Ssq[CH1 + m];	
				Phi[TIME_BLOCKS_50PC+m].real = S[CH2 + N2*m + k].real * Ssq[TIME_BLOCKS_50PC + m]; // Phi exists at each frequency bin for each channel
				Phi[TIME_BLOCKS_50PC+m].imag = S[CH2 + N2*m + k].imag * Ssq[TIME_BLOCKS_50PC + m];			
			}
			
			W_temp[0].real=0.0;
			W_temp[0].imag=0.0;
			W_temp[1].real=0.0;
			W_temp[1].imag=0.0;
			W_temp[2].real=0.0;
			W_temp[2].imag=0.0;
			W_temp[3].real=0.0;
			W_temp[3].imag=0.0;
			W_new[0].real=1.0;
			W_new[0].imag=0.0;
			W_new[1].real=0.0;
			W_new[1].imag=0.0;
			W_new[2].real=0.0;
			W_new[2].imag=0.0;
			W_new[3].real=1.0;
			W_new[3].imag=0.0;

			for(m=0;m<TIME_BLOCKS_50PC;m++) // Part of dWp(:,:,k) = (eye(nsou) - Phi*S(:,:,k)'/N)*Wp(:,:,k);
			{
				W_temp[0] = cmplx_add(W_temp[0], cmplx_mult(Phi[CH1 + m], conj(S[CH1 + N2*m + k])));
				W_temp[1] = cmplx_add(W_temp[1], cmplx_mult(Phi[CH1 + m], conj(S[CH2 + N2*m + k])));
				W_temp[2] = cmplx_add(W_temp[2], cmplx_mult(Phi[TIME_BLOCKS_50PC + m], conj(S[CH1 + N2*m + k])));
				W_temp[3] = cmplx_add(W_temp[3], cmplx_mult(Phi[TIME_BLOCKS_50PC + m], conj(S[CH2 + N2*m + k])));
			}
			
			
			W_temp[0].real = 1.0 - (W_temp[0].real / (float)TIME_BLOCKS_50PC);
			W_temp[0].imag = 0.0 - (W_temp[0].imag / (float)TIME_BLOCKS_50PC);
			W_temp[1].real = 0.0 - (W_temp[1].real / (float)TIME_BLOCKS_50PC);
			W_temp[1].imag = 0.0 - (W_temp[1].imag / (float)TIME_BLOCKS_50PC);
			W_temp[2].real = 0.0 - (W_temp[2].real / (float)TIME_BLOCKS_50PC);
			W_temp[2].imag = 0.0 - (W_temp[2].imag / (float)TIME_BLOCKS_50PC);
			W_temp[3].real = 1.0 - (W_temp[3].real / (float)TIME_BLOCKS_50PC);
			W_temp[3].imag = 0.0 - (W_temp[3].imag / (float)TIME_BLOCKS_50PC);
			
			
			// 2 by many matrix multiplied by many by 2 matrix
			W_new[0] = cmplx_add(cmplx_mult(W_temp[0], Wp[4*k + 0]), cmplx_mult(W_temp[1], Wp[4*k + 2]) );
			W_new[1] = cmplx_add(cmplx_mult(W_temp[0], Wp[4*k + 1]), cmplx_mult(W_temp[1], Wp[4*k + 3]) );
			W_new[2] = cmplx_add(cmplx_mult(W_temp[2], Wp[4*k + 0]), cmplx_mult(W_temp[3], Wp[4*k + 2]) );
			W_new[3] = cmplx_add(cmplx_mult(W_temp[2], Wp[4*k + 1]), cmplx_mult(W_temp[3], Wp[4*k + 3]) );
			
			
			/*
			for(m=0;m<NSOURCES;m++)
			{
				for(n=0;n<NSOURCES;n++)
				{
					IPhiS[m][n].real = I[m][n] - (PhiS[m][n].real/N);//I is intialised at the top
					IPhiS[m][n].imag = 0 - (PhiS[m][n].imag/N);
				}
			}*/
			// This is the "main" matrix multiplication on line 81 (original code)
		//	COMPLEX_sp_mat_mul((&IPhiS[0][0]),2,2,(&Wp[0][0][k]),2,(&dWp[0][0]));
			
			detWp = cmplx_minus(cmplx_mult(Wp[4*k + 0], Wp[4*k + 3]), cmplx_mult(Wp[4*k + 2], Wp[4*k + 3]));  // Determinate of Wp			
			dlw = dlw + (double)log(mag(detWp) + epsilon);//mag replaces the abs function
			
			// Unrolled loop - update unmixing matrix at frequency bin
			Wp[4*k + 0].real = Wp[4*k + 0].real + (mu*W_new[0].real); //mu is the learning rate
			Wp[4*k + 0].imag = Wp[4*k + 0].imag + (mu*W_new[0].imag);
			Wp[4*k + 2].real = Wp[4*k + 2].real + (mu*W_new[2].real); //mu is the learning rate
			Wp[4*k + 2].imag = Wp[4*k + 2].imag + (mu*W_new[2].imag);
			Wp[4*k + 1].real = Wp[4*k + 1].real + (mu*W_new[1].real); //mu is the learning rate
			Wp[4*k + 1].imag = Wp[4*k + 1].imag + (mu*W_new[1].imag);
			Wp[4*k + 3].real = Wp[4*k + 3].real + (mu*W_new[3].real); //mu is the learning rate
			Wp[4*k + 3].imag = Wp[4*k + 3].imag + (mu*W_new[3].imag);
		}// End frequency bin loop - line 83 in the original code
		
		
		// Sum of Ssq is worked out in the loops above before the inverse is taken.	
		//Obj = ((SumSsq/N) - dlw)/(NSOURCES*nfreq);
		Obj = (((double)SumSsq/(double)TIME_BLOCKS_50PC) - dlw)/((double)NSOURCES*(double)nfreq);
		
		
		dObj = pObj - Obj;
		pObj = Obj;
			// Update message can go here		
		comparison = fabs(dObj)/fabs(Obj);
			
		if(comparison < tol) //as Obj and dObj are real use 'abs' function
		{
			break;
		}
	} 
}

#endif /*IVA_H_*/
