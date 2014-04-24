#include "definitions.h"
#include "istft.h"

// Implements an inverse short time Fourier transform.
void istft(COMPLEX *X, float *xtime, int nfreq, int time_len, int overlap)
{
	COMPLEX buffer[1024];
	unsigned int n=0,k=0; 
	
	for(n=0;n<time_len;n+=overlap)
	{
		//reconstruct the two halfs of the FFT
		for(k=0;k<nfreq;k++)
		{
			buffer[k] = X[n + k]*hanning[k];
		}
		for(k=(nfreq-1);k>0;k--)
		{
			buffer[nfreq -1 + k] = X[n + k]*hanning[nfreq -1 + k];
		}
		
		//ifft(buffer);// This function needs to be written!
		for (k=0; k<(2*nfreq - 1); k++)
		{
			xtime[n+k] += buffer[k];
		}
	}
}