/*
 * Implements an inverse short time Fourier transform.
 * 
 * Uses TI's decimation in frequency FFT function and bit reversal. 
 * 
 * Not yet been fully debugged.
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
//#include "freqsig.h"//Test freq domain signal
#include <math.h>
#pragma DATA_ALIGN(buffer,8)
float buffer[(2*N_INT) + 8];
#pragma DATA_SECTION(buffer,".EXT_RAM")

short brev[N_INT];

void gen_w_r2(float* w, int n)
{
	int i, j=1;
 	//double pi = 4.0*atan(1.0);
 	double e = PI*2.0/n;
 	for(j=1; j < n; j <<= 1)
 	{
 		for(i=0; i < ( n>>1 ); i += j)
 		{
 			*w++ = cos((float)i*e);
 			*w++ = sin((float)i*e);
 		}
 	}
}

void istft(COMPLEX *X, float *xtime, int nfreq, int time_len, int overlap)
{
	// *X  = pointer to stft
	// *xtime =pointer to area of memory where the transform will be stored 
	// overlap = number of samples of overlap in orginal stft
	// nfreq = number of frequency bins 
	// time_len = number of samples in the original time domain signal.
	// float buffer[2*N_INT];
	unsigned int prevGIE; // Previous global interrupt flag state
	unsigned int n=0,k=0, start=4; 
	unsigned int block_ind=0;// Block index
//	complexpair *w_ptr = (complexpair*)w;// Pointer to twiddle factors

	//Not sure if the buffer array should be padded with zeros?
	for(n=0;n<8;n++)
	{
		buffer[(2*N_INT) + 8 -1 - n] = 0.0;
	}
	
	gen_w_r2(w, N_INT);
	
	bit_rev(w, N_INT>>1);/// Offending line ??? -yes!
	

	/*
	// Before the STFT a little test
	prevGIE = IRQ_globalDisable(); // Turn off global interrupts
		
 		
 	bit_rev(freq, N_INT);// BIT REVERSAL DEFINITELY GOES BEFORE!
	DSPF_sp_icfftr2_dif(freq, w, N_INT);
		
 		
	IRQ_globalRestore(prevGIE);// Restore previous gloabl interrupt state
	for(n=0; n < 2*N_INT; n++)
 	{
 		freq[n]=freq[n]/(float)N_INT;
 	}
	*/ // End of test
	// overwrite time buffer
	for(k=0; k<time_len; k++)
	{
		xtime[k]=0.0;
	}
	
	for(n=0;n<time_len;n+=overlap)// Should overlap=nfreq?
	{
		for(k=0; k<start*3; k++)//
		{
			buffer[k]=0.0;// Buffer for the buffer!
		}
		
		
		//reconstruct the two halfs of the FFT
		for(k=0; k<nfreq; k++)
		{
			buffer[2*k + 0] = X[block_ind + k].real * hamming[k]; // real
			buffer[2*k + 1] = X[block_ind + k].imag * hamming[k] * (-1); // imag
		}
		for(k=(nfreq-1); k>1; k--)// Reconstruct second half
		{
			buffer[2*(nfreq -1 + k) -2] = X[block_ind + (nfreq-1) - k].real * hamming[nfreq -2 + k]; // Real part
			buffer[2*(nfreq -1 + k) -1] = X[block_ind + (nfreq-1) - k].imag * hamming[nfreq -2 + k] * (-1); // Imag part inc. conjugate
		}
			
		block_ind += nfreq; // Go to next frequency block
		
 		prevGIE = IRQ_globalDisable(); // Turn off global interrupts
		
	
 		
 		bit_rev(buffer, N_INT);// BIT REVERSAL DEFINITELY GOES BEFORE!
		DSPF_sp_icfftr2_dif(buffer, w, N_INT); // Decimation in freq radix 2 IFFT by TI
		
	
		IRQ_globalRestore(prevGIE);// Restore previous global interrupt state
		
		
		for(k=0; k<4*(nfreq - 1); k++)
 		{
 			buffer[k]=buffer[k]/(float)N_INT;
 		}
 		
		for(k=0; k<2*(nfreq - 1); k++)
		{
			xtime[n+k] += buffer[2*k]; // Only bother with the real part
		}
	}
}
