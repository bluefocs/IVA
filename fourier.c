/*
 * Implements a short time Fourier transform and an inverse short time Fourier transform.
 * 
 * Uses TI's decimation in frequency FFT function and bit reversal for ifft. Note interrupts 
 * are turned off to use the TI DSPLIB FFT/IFFT functions.
 * 
 * Doubles the freq domain data provided to the function (i.e. recreates the other half of the 
 * STFT when creating buffer) and takes the complex conjugate.
 * 
 * Originally created 24/25-4-2014, J. Harris, Loughborough University, UK.
 * 
 * Fully debugged - 20-5-2014.
 * 
 * 
*/
#include "definitions.h"
#include "fourier.h"
#include "window.h"
#include "twiddles.h"
#include "DSPF_sp_icfftr2_dif.h"
#include "DSPF_sp_cfftr2_dit.h"
#include "DSPF_dp_cfftr2.h"
#include "csl_irq.h"
//#include "DSPF_sp_bitrev_cplx.h"
#include <math.h>

#pragma DATA_ALIGN(buffer,8)
float buffer[(2*N_INT) + 8];
#pragma DATA_SECTION(buffer,".EXT_RAM")

/*
#pragma DATA_ALIGN(buffer_dbl,16)
double buffer_dbl[(2*N_INT) + 8];
#pragma DATA_SECTION(buffer_dbl,".EXT_RAM")

#pragma DATA_ALIGN(w_dbl,16)
double w_dbl[N_INT];
#pragma DATA_SECTION(w_dbl,".EXT_RAM")
*/
float scale[CH2_T];
#pragma DATA_SECTION(scale,".EXT_RAM")

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
/*
void gen_w_r2_dbl(double* w, int n)
{
	int i, j=1;
	double pi = 4.0*atan(1.0);
	double e = pi*2.0/n;
	for(j=1; j < n; j <<= 1)
	{
		for(i=0; i < ( n>>1 ); i += j)
		{
			*w++ = cos(i*e);
			*w++ = -sin(i*e);
		}
	}
}
void bit_rev_dbl(double* x, int n)
{
	
	int i, j, k;
	double rtemp, itemp;

	j = 0;
	for(i=1; i < (n-1); i++)
	{
		k = n >> 1;
		while(k <= j)
		{
			j -= k;	
			k >>= 1;
		}
		j += k;
		if(i < j)
		{
			rtemp = x[j*2];
			x[j*2] = x[i*2];
			x[i*2] = rtemp;
			itemp = x[j*2+1];
			x[j*2+1] = x[i*2+1];
			x[i*2+1] = itemp;
		}
	}
}
*/

void istft(COMPLEX *X, float *xtime, int nfreq, int time_len, int overlap)
{
	// *X  = pointer to stft
	// *xtime =pointer to area of memory where the inverse transform will be stored 
	// overlap = number of samples of overlap in orginal stft
	// nfreq = number of frequency bins 
	// time_len = number of samples in the original time domain signal.

//	unsigned int prevGIE; // Previous global interrupt flag state
	unsigned int n=0,k=0, start=4; 
	unsigned int block_ind=0, step = N_INT;// Block index
	const float fftlen_inv = 1/(float)N_INT;
//	complexpair *w_ptr = (complexpair*)w;// Pointer to twiddle factors


	step = step - overlap;

	// Variables for scaling
	//float scale=0.0, w0=0.0;

	//Not sure if the buffer array should be padded with zeros?
	for(n=0;n<8;n++)
	{
		buffer[(2*N_INT) + 8 -1 - n] = 0.0;
	}
	
	gen_w_r2(w, N_INT);
	
	bit_rev(w, N_INT>>1);/// Offending line ??? -yes!
	

	/*
	// Before the STFT a little test - left commented out incase one needs a quick sanity check ...
	prevGIE = IRQ_globalDisable(); // Turn off global interrupts
		
 		
 	bit_rev(freq, N_INT);// BIT REVERSAL DEFINITELY GOES BEFORE!
	DSPF_sp_icfftr2_dif(freq, w, N_INT);
		
 		
	IRQ_globalRestore(prevGIE);// Restore previous gloabl interrupt state
	for(n=0; n < 2*N_INT; n++)
 	{
 		freq[n]=freq[n]*fftlen_inv;
 	}
	*/ // End of test
	// overwrite time buffer
	for(k=0; k<time_len; k++)
	{
		xtime[k]=0.0;
		scale[k]=0.0;
	}
	
	for(n=0; n<(time_len-(nfreq-1)*2); n+=step)// Should overlap=nfreq?
	{
		for(k=0; k<start*3; k++)//
		{
			buffer[k]=0.0;// Buffer for the buffer!
		}
		
		
		//reconstruct the two halfs of the FFT
		for(k=0; k<nfreq; k++)
		{
			buffer[2*k + 0] = X[block_ind + k].real; // real
			buffer[2*k + 1] = X[block_ind + k].imag; // imag - take the conj here as this accounts for not taking at the beginning of the funciton
		}
		for(k=0; k<(nfreq-2); k++)// Reconstruct second half, (nfreq-2) starting point - we don't want the first val (real) of the fft again
		{
			buffer[2 * (nfreq + k)] = X[block_ind + nfreq - 2 - k].real; // Real part
			buffer[2*(nfreq + k)+1] = X[block_ind + nfreq - 2 - k].imag * (-1); // Imag part inc. conjugate
		}
			
		block_ind += nfreq; // Go to next time block
		
 		
		
	
 		bit_rev(buffer, N_INT);// BIT REVERSAL DEFINITELY GOES BEFORE! (IFFT)
 		//prevGIE = IRQ_globalDisable(); // Turn off global interrupts
		DSPF_sp_icfftr2_dif_c(buffer, w, N_INT); // Decimation in freq radix 2 IFFT by TI	
		//IRQ_globalRestore(prevGIE);// Restore previous global interrupt state
		
		

 		
		for(k=0; k<2*(nfreq - 1); k++) // Overlap add method
		{
			xtime[n+k] += buffer[2*k] * hanning[k] * fftlen_inv; // Only bother with the real part
			
			scale[n+k] += hanning[k] * hanning[k]; // Work out weird scally factor thing
		}			
	}
	
	
	
	for(k=0; k<2*(nfreq - 1); k++) // Overlap add method
	{
		xtime[n+k] += buffer[2*k] * hanning[k]; // Only bother with the real part		
		scale[n+k] += hanning[k] * hanning[k]; // Work out weird scally factor thing
	}
	
	
	
/*	for(k=0; k<2*(nfreq - 1); k++) // Work out last little bit of the scally thing
	{
		scale[n+k] += hamming[k] * hamming[k]; // Only bother with the real part
	}*/ 
	// This scaling is included in other implementations of istft
	/*for(k=0; k<2*(nfreq - 1); k++) // Calculate scaling factor
	{
		w0 += hamming[k]*hamming[k];
	}

	scale = overlap/w0;*/
	scale[0]=scale[1];// Ensures that the first value isn't 0 <- bit of a bodge
	for(k=0; k<time_len; k++) // Calculate scaling factor
	{
		xtime[k] = xtime[k] / scale[k];
	}
}


void stft(COMPLEX *X, float *xtime, int nfreq, int time_len, int overlap)
{
	/* 
	*
	* Uses hanning window, divided by length then multiplied by 2.
	*/
	unsigned int prevGIE; // Previous global interrupt flag state
	unsigned int n=0, m=0, step=0, overlap_pos=0;
	unsigned short k=0;
	
	//unsigned int no_time_frames <--- was going to be used to introduce zero padding
	step = 2*(nfreq-1) - overlap;
	
	for(n=0;n<8;n++)
	{
		buffer[(2*N_INT) + 8 -1 - n] = 0.0;
	}
	// Set up twiddle factors to be sure that they work with the TI function
	//gen_w_r2_dbl(w_dbl, N_INT);
	gen_w_r2(w, N_INT);
	bit_rev(w, N_INT>>1);/// SINGLE (FLOAT) Offending line ??? -yes!

	

	// This loop effectively pads the time domain signal with zeros
	m=0;//Important!
	for(n=1; n<((nfreq-1)*2/step); n++) // N/2 for 50% overlapping 
	{	
		overlap_pos = step*((nfreq-1)*2/step-n);
		
		for(k = 0; k < overlap_pos; k++)
		{
			buffer[2 * k] = 0.0;//hanning[k] * xtime[n+k];
			buffer[2*k+1] = 0.0;
		}
		for(k = overlap_pos; k<(nfreq-1)*2; k++)
		{
			buffer[2 * k] = hanning[k] * xtime[k - overlap_pos];
			buffer[2*k+1] = 0.0;
		}
		// BUT! memcpy seems to be more efficient - this way you can't use the window
		//memcpy(&buffer, &X[n], N*sizeof(complexpair)); // Copy full 1024 time domain points
		
		// Calculate 1024 point FFT on current buffers
		prevGIE = IRQ_globalDisable(); // Turn off global interrupts
		
		DSPF_sp_cfftr2_dit_c(buffer, w, N_INT);
		bit_rev(buffer, N_INT);
		
		//DSPF_sp_cfftr2_dit(buffer, w, N_INT); // Pre-compiled version of above function 
		//DSPF_dp_cfftr2(N_INT, buffer_dbl, w_dbl, 1);// Double precision experiment - didn't work!
		//bit_rev_dbl(buffer_dbl, N_INT); // Bit reversal goes afterwards for the forward fft
		IRQ_globalRestore(prevGIE);
		
		memcpy(&X[m], &buffer[0], N*sizeof(complexpair)); // N is half +1 of N_INT
		m += nfreq;
	}




	
	//m=0;//Important!
	//for(n=0; n<((N_INT*TIME_BLOCKS_INT)-(N_INT/2)); n+=step) // N/2 for 50% overlapping 
	for(n=0; n<((N_INT*TIME_BLOCKS_INT)-N_INT); n+=step) // N/2 for 50% overlapping 
	{	
		
		// In order to implement the window you need to loop around every value and multiply it by the relevant coefficient 
		for(k=0;k<N_INT;k++)
		{
			buffer[2 * k] = hanning[k] * xtime[n+k];
			buffer[2*k+1]= 0.0;
		//	buffer_dbl[2 * k] = (double)hanning[k] * (double)xtime[n+k];
		//	buffer_dbl[2*k+1]= 0.0;
		}
		// BUT! memcpy seems to be more efficient - this way you can't use the window
		//memcpy(&buffer, &X[n], N*sizeof(complexpair)); // Copy full 1024 time domain points
		
		// Calculate 1024 point FFT on current buffers
		prevGIE = IRQ_globalDisable(); // Turn off global interrupts
		
		DSPF_sp_cfftr2_dit_c(buffer, w, N_INT);
		bit_rev(buffer, N_INT);
		
		//DSPF_sp_cfftr2_dit(buffer, w, N_INT); // Pre-compiled version of above function 
		//DSPF_dp_cfftr2(N_INT, buffer_dbl, w_dbl, 1);// Double precision experiment - didn't work!
		//bit_rev_dbl(buffer_dbl, N_INT); // Bit reversal goes afterwards for the forward fft
		IRQ_globalRestore(prevGIE);
		
		/*
		for(k=0;k<N_INT;k++)// These lines were for the double test
		{
			buffer[2 * k] = (float)buffer_dbl[2 * k];
			buffer[2*k+1] = (float)buffer_dbl[2*k+1];
		}*/
		
		memcpy(&X[m], &buffer[0], N*sizeof(complexpair)); // N is half +1 of N_INT
		m += nfreq;
	}
	
	//Overlap with 000s at the end
	for(n=0; n<((nfreq-1)*2/step); n++) // N/2 for 50% overlapping 
	{	
		overlap_pos = step*((nfreq-1)*2/step-n);
		
		for(k = 0; k < step*(n+1); k++)
		{
			buffer[2 * k] = hanning[k] * xtime[CH2_T - overlap + n*step + k];
			buffer[2*k+1] = 0.0;
		}
		for(k = step*(n+1); k<(nfreq-1)*2; k++)
		{
			buffer[2 * k] = 0.0;//hanning[k] * xtime[n+k];
			buffer[2*k+1] = 0.0;
		}
		// BUT! memcpy seems to be more efficient - this way you can't use the window
		//memcpy(&buffer, &X[n], N*sizeof(complexpair)); // Copy full 1024 time domain points
		
		// Calculate 1024 point FFT on current buffers
		prevGIE = IRQ_globalDisable(); // Turn off global interrupts
		
		DSPF_sp_cfftr2_dit_c(buffer, w, N_INT);
		bit_rev(buffer, N_INT);
		
		//DSPF_sp_cfftr2_dit(buffer, w, N_INT); // Pre-compiled version of above function 
		//DSPF_dp_cfftr2(N_INT, buffer_dbl, w_dbl, 1);// Double precision experiment - didn't work!
		//bit_rev_dbl(buffer_dbl, N_INT); // Bit reversal goes afterwards for the forward fft
		IRQ_globalRestore(prevGIE);
		
		memcpy(&X[m], &buffer[0], N*sizeof(complexpair)); // N is half +1 of N_INT
		m += nfreq;
	}	
}

void DSPF_sp_cfftr2_dit_c(float* x, float* w, short n)
{                                                                   
	short n2, ie, ia, i, j, k, m;                                    
	float rtemp, itemp, c, s;                                        
                                                                           
	n2 = n;                                                          
	ie = 1;                                                          
                                                                           
	for(k=n; k > 1; k >>= 1)                                         
	{                                                                
		n2 >>= 1;                                                     
		ia = 0;                                                       
		for(j=0; j < ie; j++)                                         
		{                                                             
			c = w[2*j];                                                
			s = w[2*j+1];                                              
			for(i=0; i < n2; i++)                                      
			{                                                          
				m = ia + n2;                                            
				rtemp     = c * x[2*m]   + s * x[2*m+1];                
				itemp     = c * x[2*m+1] - s * x[2*m];                  
				x[2*m]    = x[2*ia]   - rtemp;                          
				x[2*m+1]  = x[2*ia+1] - itemp;                          
				x[2*ia]   = x[2*ia]   + rtemp;                          
				x[2*ia+1] = x[2*ia+1] + itemp;                          
				ia++;                                                   
			}                                                          
			ia += n2;                                                  
		}                                                             
		ie <<= 1;                                                     
	}                                                                
}  
void DSPF_sp_icfftr2_dif_c(float    * x, float * w, short n)
{

/* ======================================================================== */
/* ======================================================================== */
     short n2, ie, ia, i, j, k, m;
     float rtemp, itemp, c, s;

     n2 = 1;
     ie = n;
     for(k=n; k > 1; k >>= 1)
     {
        ie >>= 1;
        ia = 0;
        for(j=0; j < ie; j++)
        {
           c = w[2*j];
           s = w[2*j+1];
           for(i=0; i < n2; i++)
           {
              m = ia + n2;
              rtemp     = x[2*ia]   - x[2*m];
              x[2*ia]   = x[2*ia]   + x[2*m];
              itemp     = x[2*ia+1] - x[2*m+1];
              x[2*ia+1] = x[2*ia+1] + x[2*m+1];
              x[2*m]    = c*rtemp   - s*itemp;
              x[2*m+1]  = c*itemp   + s*rtemp;
              ia++;
           }
           ia += n2;
        }
        n2 <<= 1;
      }
}
