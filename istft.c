/*
 * Implements an inverse short time Fourier transform.
 * 
 * Uses TI's decimation in frequency FFT function and bit reversal. 
 * 
 * Not let been fully debugged.
 * 
 * Originally created 24/25-4-2014, J. Harris, Loughborough University, UK
 * 
 * 
*/
#include "definitions.h"
#include "istft.h"
#include "hamming.h"
#include "twiddles.h"
#include "DSPF_sp_icfftr2_dif.h"
#include "csl_irq.h"
#include "DSPF_sp_bitrev_cplx.h"

void istft(COMPLEX *X, float *xtime, int nfreq, int time_len, int overlap)
{
	// *X  = pointer to stft
	// *xtime =pointer to area of memory where the transform will be stored 
	// overlap = number of samples of overlap in orginal stft
	// nfreq = number of frequency bins 
	// time_len = number of samples in the original time domain signal.
	float buffer[2*N_INT];
	unsigned int prevGIE; // Previous global interrupt flag state
	unsigned int n=0,k=0; 
	short brev[N_INT];
	//float xtest[16] = {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
	//short brev[8] = {0,0,0,0,0,0,0,0};

	for(n=0; n<N_INT; n++)
	{
		brev[n] = 0;
	}
	
	bitrev_index(brev,N_INT);
	//DSPF_sp_bitrev_cplx((double*)xtest, brev, N_INT);//N_INT=1024
	
	for(n=0;n<time_len;n+=overlap)
	{
		//reconstruct the two halfs of the FFT
		for(k=0; k<nfreq; k++)
		{
			buffer[2*k + 0] = X[n + k].real*hamming[k]; // real
			buffer[2*k + 1] = X[n + k].imag*hamming[k]; // imag
		}
		for(k=(nfreq-1); k>1; k--)// Reconstruct second half
		{
			buffer[(nfreq -1)*2 + 2*k -2] = X[n + (nfreq-1) - k].real * hamming[nfreq -2 + k]; // Need to take the conjugate here
			buffer[(nfreq -1)*2 + 2*k -1] = X[n + (nfreq-1) - k].imag * hamming[nfreq -2 + k];
		}
		


		
 		prevGIE = IRQ_globalDisable(); // Turn off global interrupts
		
		// Bit reverse  frequency domain data
		DSPF_sp_bitrev_cplx((double*)buffer, brev, N_INT);//N_INT=1024
		DSPF_sp_icfftr2_dif(buffer, w, ((nfreq*2)-1));
 
		IRQ_globalRestore(prevGIE);// Restore previous gloabl interrupt state
		for (k=0; k<2*(nfreq - 1); k++)
		{
			xtime[n+k] += buffer[2*k]; // Only bother with the real part
		}
	}
}
