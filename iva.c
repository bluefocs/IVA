#include "iva.h"
#include "definitions.h"
#include "additional_math.h"
#include "DSPF_sp_mat_mul_cplx.h"

COMPLEX S[N * NSOURCES * TIME_BLOCKS]; //
#pragma DATA_SECTION(S,".EXT_RAM")


void iva(COMPLEX *Xp, COMPLEX *Wp, unsigned short nfreq)
{
	const float recip_TIME_BLOCKS = 1.0 / (float)TIME_BLOCKS;
	unsigned short k=0, m=0, mat_k=0, iter=0;
	unsigned int ch1_i=0, ch2_i=0;
	const unsigned short maxiter=600;
	const float mu=0.1;
//	COMPLEX detWp;
	float Ssq[TIME_BLOCKS * NSOURCES];
//	double epsilon = 0.00001;
	COMPLEX W_temp[4], W_new[4];//, detWp;	
	COMPLEX Phi[TIME_BLOCKS * NSOURCES];
//	double SumSsq=0.0, dObj=0.0, pObj=0.0, Obj=0.0, dlw=0.0, tol = 0.000001, comparison=0.0;
	
	// Initialise Ssq
	for  (m=0; m<(TIME_BLOCKS * NSOURCES); m++)
	{
		Ssq[m]=0.0;
	}
	
	#pragma MUST_ITERATE(600,600)
	for(iter=0;iter<maxiter;iter++)
	{
		//dlw = 0;// Used for cost function value		
		#pragma MUST_ITERATE(513,513)
		for(k=0;k<N;k++)
		{
			#pragma MUST_ITERATE(TIME_BLOCKS,TIME_BLOCKS)
			for(m=0; m<TIME_BLOCKS; m++)// 2 by many matrix multiplied by many by 2 matrix
			{
				ch1_i = N*m + k; // Channel 1 index
				ch2_i = CH2 + ch1_i; // Channel 2 index
				mat_k = 4*k;// matrix index
				cmplx_mult_add(Wp[mat_k + 0], Xp[ch1_i], Wp[mat_k + 1], Xp[ch2_i], &S[ch1_i].real, &S[ch1_i].imag);
				cmplx_mult_add(Wp[mat_k + 2], Xp[ch1_i], Wp[mat_k + 3], Xp[ch2_i], &S[ch2_i].real, &S[ch2_i].imag);
				//S[CH1 + i] = cmplx_add(cmplx_mult(Wp[4*k + 0], Xp[CH1 + i]), cmplx_mult(Wp[4*k + 1], Xp[CH2 + i]));
				//S[CH2 + i] = cmplx_add(cmplx_mult(Wp[4*k + 2], Xp[CH1 + i]), cmplx_mult(Wp[4*k + 3], Xp[CH2 + i]));
			}
		}
		
		
		
//		SumSsq=0.0;// Used for cost function value
		// Calculate score function function - derived from the multivariate Gaussian distribution function.
		for(m=0;m<TIME_BLOCKS;m++) // Summnation loop - Can this be sped up ? and be done with double precision?
		{
			for(k=0;k<nfreq;k++)
			{
				ch1_i = N*m + k;		// Channel 1 index
				ch2_i = CH2 + ch1_i; 	// Channel 2 index		
				Ssq[ m ] += mag_sqrd(S[ch1_i]);	
				Ssq[TIME_BLOCKS+m] += mag_sqrd(S[ch2_i]);	
				//Ssq[ m ] += ((mag(S[ch1_i])) * (mag(S[ch1_i])));			//pow(mag(S[CH1 + i]), 2.0);
				//Ssq[TIME_BLOCKS+m] += ((mag(S[CH2 + i])) * (mag(S[CH2 + i])));//pow(mag(S[CH2 + i]), 2.0);

				// Non - assembly way of doing it 
				//Ssq[ m ] += ((mag(S[CH1 + i])) * (mag(S[CH1 + i])));			//pow(mag(S[CH1 + i]), 2.0);
				//Ssq[TIME_BLOCKS+m] += ((mag(S[CH2 + i])) * (mag(S[CH2 + i])));//pow(mag(S[CH2 + i]), 2.0);
				// Use TI's optimised fastmath library
				//Ssq[ m ] += powsp(mag(S[CH1 + N2*m + k]),2.0);
				//Ssq[TIME_BLOCKS_50PC+m] += powsp(mag(S[CH2 + N2*m + k]), 2.0);
			}
			
			
			Ssq[ m ] = FastInvSqrt(Ssq[ m ]); // In the future change ^0.5 to ^0.666. Important line! 
			Ssq[TIME_BLOCKS+m] = FastInvSqrt(Ssq[TIME_BLOCKS+m]); // Channel 2
			
			
			/*
			// The 5 lines below can be sped up by using the fast inverse sqrt function
			Ssq[ m ] = sqrt(Ssq[ m ]); // In the future change ^0.5 to ^0.666. Important line! 
			Ssq[TIME_BLOCKS+m] = sqrt(Ssq[TIME_BLOCKS+m]); // Channel 2
			
			// Calculate the sum of all the values in Ssq before the inverse is taken, this is used in the break condition below
			SumSsq += Ssq[ m ] + Ssq[TIME_BLOCKS+m];
			
			// This 'inversion comes from the derivative of the cost function, G'(ri)/ri (See Yanfeng's EUSIPCO paper)
			// Does the inversion work in this loop? 
			Ssq[ m ] = 1.0/(Ssq[ m ] + epsilon);	// Channel 1 - Ssq1 in MATLAB code
			Ssq[TIME_BLOCKS+m]=1.0/(Ssq[TIME_BLOCKS+m] + epsilon);	//Channel 2
			*/
		}
		

		/*
		// This 'inversion comes from the derivative of the cost function, G'(ri)/ri (See Yanfeng's EUSIPCO paper)
		for(m=0;m<TIME_BLOCKS;m++) // Take the summnation of 
		{
			Ssq[ m ] = 1.0/(Ssq[ m ] + epsilon);	// Channel 1 - Ssq1 in MATLAB code
			Ssq[TIME_BLOCKS+m]=1.0/(Ssq[TIME_BLOCKS+m] + epsilon);	//Channel 2
		}
		*/

		

		for(k=0;k<nfreq;k++)
		{
			mat_k = 4*k;// Matrix index for unmixing matrices, etc.
			//Calculate multivariate score function and gradients
			for(m=0; m<TIME_BLOCKS; m++)
			{
				// Unrolled the source loop
				ch1_i = N*m + k; // Channel 1 index
				ch2_i = ch1_i + CH2; // Channel 2 index
				Phi[ CH1 + m ].real = S[ch1_i].real * Ssq[CH1 + m]; // Phi exists at each frequency bin for each channel
				Phi[ CH1 + m ].imag = S[ch1_i].imag * Ssq[CH1 + m];	
				Phi[TIME_BLOCKS+m].real = S[ch2_i].real * Ssq[TIME_BLOCKS + m]; // Phi exists at each frequency bin for each channel
				Phi[TIME_BLOCKS+m].imag = S[ch2_i].imag * Ssq[TIME_BLOCKS + m];			
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

			for(m=0;m<TIME_BLOCKS;m++) // Part of dWp(:,:,k) = (eye(nsou) - Phi*S(:,:,k)'/N)*Wp(:,:,k);
			{
				ch1_i = N*m + k; // Channel 1 index
				ch2_i = ch1_i + CH2; // Channel 2 index
				
				W_temp[0] = cmplx_add(W_temp[0], cmplx_mult(Phi[CH1 + m], conj(S[ch1_i])));
				W_temp[1] = cmplx_add(W_temp[1], cmplx_mult(Phi[CH1 + m], conj(S[ch2_i])));
				W_temp[2] = cmplx_add(W_temp[2], cmplx_mult(Phi[TIME_BLOCKS + m], conj(S[ch1_i])));
				W_temp[3] = cmplx_add(W_temp[3], cmplx_mult(Phi[TIME_BLOCKS + m], conj(S[ch2_i])));
			}
			
			
			W_temp[0].real = 1.0 - (W_temp[0].real * recip_TIME_BLOCKS);// / (float)TIME_BLOCKS);
			W_temp[0].imag = 0.0 - (W_temp[0].imag * recip_TIME_BLOCKS);// / (float)TIME_BLOCKS);
			W_temp[1].real = 0.0 - (W_temp[1].real * recip_TIME_BLOCKS);// / (float)TIME_BLOCKS);
			W_temp[1].imag = 0.0 - (W_temp[1].imag * recip_TIME_BLOCKS);// / (float)TIME_BLOCKS);
			W_temp[2].real = 0.0 - (W_temp[2].real * recip_TIME_BLOCKS);// / (float)TIME_BLOCKS);
			W_temp[2].imag = 0.0 - (W_temp[2].imag * recip_TIME_BLOCKS);// / (float)TIME_BLOCKS);
			W_temp[3].real = 1.0 - (W_temp[3].real * recip_TIME_BLOCKS);// / (float)TIME_BLOCKS);
			W_temp[3].imag = 0.0 - (W_temp[3].imag * recip_TIME_BLOCKS);// / (float)TIME_BLOCKS);
			
			// 2 2 by many matrix multiplied by many by 2 matrix
			cmplx_mult_add(W_temp[0], Wp[mat_k + 0], W_temp[1], Wp[mat_k + 2], &W_new[0].real, &W_new[0].imag);
			cmplx_mult_add(W_temp[0], Wp[mat_k + 1], W_temp[1], Wp[mat_k + 3], &W_new[1].real, &W_new[1].imag);
			cmplx_mult_add(W_temp[2], Wp[mat_k + 0] ,W_temp[3], Wp[mat_k + 2], &W_new[2].real, &W_new[2].imag);
			cmplx_mult_add(W_temp[2], Wp[mat_k + 1], W_temp[3], Wp[mat_k + 3], &W_new[3].real, &W_new[3].imag);
			// 2 by many matrix multiplied by many by 2 matrix
			//W_new[0] = cmplx_add(cmplx_mult(W_temp[0], Wp[4*k + 0]), cmplx_mult(W_temp[1], Wp[4*k + 2]) );
			//W_new[1] = cmplx_add(cmplx_mult(W_temp[0], Wp[4*k + 1]), cmplx_mult(W_temp[1], Wp[4*k + 3]) );
			//W_new[2] = cmplx_add(cmplx_mult(W_temp[2], Wp[4*k + 0]), cmplx_mult(W_temp[3], Wp[4*k + 2]) );
			//W_new[3] = cmplx_add(cmplx_mult(W_temp[2], Wp[4*k + 1]), cmplx_mult(W_temp[3], Wp[4*k + 3]) );
			
			
			
		//	detWp = cmplx_minus(cmplx_mult(Wp[4*k + 0], Wp[4*k + 3]), cmplx_mult(Wp[4*k + 2], Wp[4*k + 3]));  // Determinate of Wp			
		//	dlw = dlw + (double)log(mag(detWp) + epsilon);//mag replaces the abs function
			
			// Unrolled loop - update unmixing matrix at frequency bin
			Wp[mat_k + 0].real = Wp[mat_k + 0].real + (mu*W_new[0].real); //mu is the learning rate
			Wp[mat_k + 0].imag = Wp[mat_k + 0].imag + (mu*W_new[0].imag);
			Wp[mat_k + 2].real = Wp[mat_k + 2].real + (mu*W_new[2].real); //mu is the learning rate
			Wp[mat_k + 2].imag = Wp[mat_k + 2].imag + (mu*W_new[2].imag);
			Wp[mat_k + 1].real = Wp[mat_k + 1].real + (mu*W_new[1].real); //mu is the learning rate
			Wp[mat_k + 1].imag = Wp[mat_k + 1].imag + (mu*W_new[1].imag);
			Wp[mat_k + 3].real = Wp[mat_k + 3].real + (mu*W_new[3].real); //mu is the learning rate
			Wp[mat_k + 3].imag = Wp[mat_k + 3].imag + (mu*W_new[3].imag);
		}// End frequency bin loop - line 83 in the original code
		
		
		// Sum of Ssq is worked out in the loops above before the inverse is taken.	
		
		
		/*
		Obj = (((double)SumSsq/(double)TIME_BLOCKS) - dlw)/((double)NSOURCES*(double)nfreq);
		
		
		dObj = pObj - Obj; // Work out change in Obj
		pObj = Obj; // Current Obj becomes 'previous' Obj
			// Update message can go here		
		comparison = fabs(dObj)/fabs(Obj);//as Obj and dObj are real use 'fabs' function
			
		if(comparison < tol) // Is the change in Obj divided by the current Obj lower than a tolerance?
		{
			break;
		}*/
	} 
}
