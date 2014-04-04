#ifndef IVA_H_
#define IVA_H_
#include "definitions.h"

COMPLEX S[N2 * 2 * TIME_BLOCKS_50PC]; //

void iva(COMPLEX *Xp, COMPLEX *Wp, unsigned short nfreq)
{
	short dlw=0,k=0, m=0;
	unsigned short maxiter=1000, iter=0;
	float mu=0.1;
	COMPLEX detWp;
	float Ssq[TIME_BLOCKS_50PC * NSOURCES];
	float epsilon = 0.000001;
	COMPLEX W_temp[4], W_new[4];//, detWp;	
	COMPLEX Phi[TIME_BLOCKS_50PC * NSOURCES];
	float SumSsq, dObj, pObj, Obj;
	float tol=1e-6;
	
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
		
		// Calculate non-linear function?
		for(m=0;m<TIME_BLOCKS_50PC;m++) // Summnation loop
		{
			for(k=0;k<nfreq;k++)
			{
				Ssq[CH1 + m] += pow(mag(S[CH1 + N2*m + k]), 2.0);
				Ssq[CH2 + m] += pow(mag(S[CH2 + N2*m + k]), 2.0);
			}
		}
		for(m=0;m<TIME_BLOCKS_50PC;m++) // Take the sumnation of 
		{
			Ssq[ m ] = 1.0/(Ssq[CH1 + m] + epsilon);	// Channel 1 
			Ssq[TIME_BLOCKS_50PC + m]=1.0/(Ssq[CH2 + m] + epsilon);	//Channel 2
		}
		
		

		for(k=0;k<nfreq;k++)
		{
			//Calculate multivariate score function and gradients
			for(m=0; m<TIME_BLOCKS_50PC; m++)
			{
				//for(n=0;n<NSOURCES;n++)
				// Unrolled the source loop
					
				Phi[CH1 + m].real = S[CH1 + N2*m + k].real * Ssq[CH1 + m]; // Phi exists at each frequency bin for each channel
				Phi[CH1 + m].imag = S[CH1 + N2*m + k].imag * Ssq[CH1 + m];	
				Phi[CH2 + m].real = S[CH2 + N2*m + k].real * Ssq[CH2 + m]; // Phi exists at each frequency bin for each channel
				Phi[CH2 + m].imag = S[CH2 + N2*m + k].imag * Ssq[CH2 + m];			
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
				W_temp[0] = cmplx_add(W_temp[0], cmplx_mult(Phi[CH1 + m], S[CH2 + N2*m + k]));
				W_temp[1] = cmplx_add(W_temp[1], cmplx_mult(Phi[CH1 + m], S[CH1 + N2*m + k]));
				W_temp[2] = cmplx_add(W_temp[2], cmplx_mult(Phi[CH2 + m], S[CH1 + N2*m + k]));
				W_temp[3] = cmplx_add(W_temp[3], cmplx_mult(Phi[CH2 + m], S[CH2 + N2*m + k]));
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
			dlw = dlw + log(mag(detWp) + epsilon);//mag replaces the abs function
			
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
		
		
		SumSsq = Ssq[0] + Ssq[1] + Ssq[2] + Ssq[3];// Bit of a rubbish way of doing the summnation but for now it's easier to read
		//SumSsq.imag = Ssq[0][0].imag + Ssq[0][1].imag + Ssq[1][0].imag + Ssq[1][1].imag;
			
			
		Obj = ((SumSsq/N) - dlw)/(NSOURCES*nfreq);
		
		dObj = pObj - Obj;
		pObj = Obj;
			// Update message can go here		
			
			
		if(abs(dObj)/abs(Obj)<tol) //as Obj and dObj are real use 'abs' function
		{
			break;
		}
	} 
}

#endif /*IVA_H_*/
