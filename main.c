/*
 * main.c - written by Jack Harris Feburary-March-April 2014, Loughborough, UK.
 * 
 * 
 * 
 * LED 1) Reads in approx 5 seconds of audio data in progress.
 * LED 2) Indicates 1024 point FFT / STFT is complete with 50% overlapping with hanning 
 * 			window (see note).
 * LED 3) PCA complete.
 * LED 4) IVA complete. Though feedthrough audio continues via the interrupt.
 * 
 * Interrupt - continuous 'feedthrough' audio system continuously sampling audio data
 * from the line in and outputting it to the headphones.
 * git st
 * Twiddle values for the FFT are generated in MATLAB (generate_twiddles.m) due to 
 * previous problems with the cos (cosine) function with the TI toolchain.
 * 
 * The STFT is of length 1024 and uses a hamming window (the values have been precalculated
 * and placed in hamming.h). The memory works by placing the time domain signal in memory 
 * space of 40*1024 (as if there were going to be 40 non-overlapping frequency bins), but 
 * places the result of the STFT in memory space 79*1024, the increase in blocks is due to 
 * the 50% overlap. 
 * 
 * The PCA and whitening has been fully debugged. Returns identity matrix when two 
 * corresponding frequency bins are compared with the cov() function in MATLAB
 * (i.e. they are uncorrelated).
 * 
 * IVA has been been debugged. This has been verified with the output of the IVA 
 * algorithm in MATLAB (provided by Taesu Kim). The exit clause needs to be investigated,
 * this version runs too long/too short. However returns the correct Wp (or near enough).
 * 
 * The program can suffer from NaN propogating through each stage, this is POSSIBLY 
 * caused by the FFT code from Chassaing's book and dealing with large input values, 
 * in future versions scaling needs to be considered. Possibly turn down the volume 
 * of whatever is driving the audio in. This has been potential solved on 14-4-2014 by 
 * changing the outdata variable from unsigned to signed (!)
 * 
 * Losely based on the original file in Rulph Chassaing's TMS320C6713 DSK book:
 * frames.c - basic illustration of triple-buffered
 *  N sample frame-based processing
 * 
 * TO DO LIST:
 *  - Inverse STFT still needs to be implemented. - completed 25-4-2014 (further verification may be required)
 *  - (Possibly) expand frequency bins from 512 to 1024 (don't throw away half the 
 * 		frequency bins to speed up processing) - completed 14-4-2014
 * 	- Implement a complex multiply function in assembley. completed - 26-4-2014.
 * 	- Multiply the unmixing matrices by the estimated sources at each frequency bin.
 *  - Need to rewrite the exit clause for the IVA algorithm (checking the cost function).
 *  - Write FastInvSqrt function in linear assembly. 
 */
#include "dsk6713.h"
#include "DSK6713_AIC23.h"	//codec-DSK interface support
#include "fft.h"
#include "additional_math.h"
#include "iva.h"
#include "definitions.h"
#include "twiddles.h"
#include "hamming.h"
#include "istft.h"
#include <csl_cache.h>
//#include <utility.h>
// Inlcudes for ifft test
#include "DSPF_sp_icfftr2_dif.h"
#include "DSPF_sp_bitrev_cplx.h"

Uint32 fs=DSK6713_AIC23_FREQ_8KHZ;	//set sampling rate
DSK6713_AIC23_CodecHandle hCodec; // codec handle declaration
Uint16 inputsource=DSK6713_AIC23_INPUT_LINE; // select input source


union complexdata X[NSOURCES*STFT_SIZE];// Used to store time data, freq domain data and whitened data
float x[NSOURCES * TIME_BLOCKS_INT * N_INT];
float x_rec[TIME_BLOCKS_INT * N_INT];
unsigned short buffercount = 0;            //number of new input samples in buffer 
static unsigned short t = 0;
complexpair *X_ptr[TIME_BLOCKS];
COMPLEX *Xstart_ptr;


// Buffers for the FFT
complexpair buffer1[N_INT],buffer2[N_INT];	// Buffers for the FFTs of length 1024

// IVA variables - placed here so that the code works
COMPLEX Q[N*4];// 2*2 matrix at each freq bin
COMPLEX Wp[N*4];
#pragma DATA_SECTION(X,".EXT_RAM")
#pragma DATA_SECTION(Q,".EXT_RAM")
#pragma DATA_SECTION(Wp,".EXT_RAM")
#pragma DATA_SECTION(buffer1,".EXT_RAM")
#pragma DATA_SECTION(buffer2,".EXT_RAM")
#pragma DATA_SECTION(x,".EXT_RAM")
#pragma DATA_SECTION(x_rec,".EXT_RAM")
//#pragma DATA_SECTION(DSPF_sp_icfftr2_dif,".EXT_RAM")

// Function prototypes for any assembly routines
//extern COMPLEX cmplx_mult_sp(COMPLEX x, COMPLEX y, float *real, float *imag);

void tw_genSPxSPfft(float * w, int n)                                      
{                                                                  
     int i, j, k;                                                        
     double x_t, y_t, theta1, theta2, theta3;                         
                                                                       
     for (j=1, k=0; j <= n>>2; j = j<<2)                              
     {                                                                
         for (i=0; i < n>>2; i+=j)                                    
         {                                                            
             theta1 = 2*PI*i/n;                                       
             x_t = cos(theta1);                                       
             y_t = sin(theta1);                                       
             w[k]   =  (float)x_t;                                    
             w[k+1] =  (float)y_t;                                    
                                                                      
             theta2 = 4*PI*i/n;                                       
             x_t = cos(theta2);                                       
             y_t = sin(theta2);                                       
             w[k+2] =  (float)x_t;                                    
             w[k+3] =  (float)y_t;                                    
                                                                      
             theta3 = 6*PI*i/n;                                       
             x_t = cos(theta3);                                       
             y_t = sin(theta3);                                       
             w[k+4] =  (float)x_t;                                    
             w[k+5] =  (float)y_t;                                    
             k+=6;                                                    
         }                                                            
     }                                                                
}


interrupt void c_int11(void)      //ISR
{
  	union {Uint32 uint; short channel[2];} outdata;
  	
	DSK6713_LED_on(0);		
	
	outdata.uint = input_sample(); // 
	output_sample(outdata.uint);
	
	//Original way of reading in data
	//(input_ptr + buffercount)->numbers[IMAG] = 0.0; 
	//(input_ptr + buffercount++)->numbers[REAL] = (float)(outdata.channel[LEFT]); 
		
	if(t<(TIME_BLOCKS_INT))
	{	
		//(X_ptr[t] + buffercount)->numbers[IMAG] = 0.0; 
		//(X_ptr[t] + buffercount++)->numbers[REAL] = (float)(outdata.channel[LEFT]);	
		
		//X[CH1 + (t*N) + buffercount].numbers[REAL] = (float)(input_left_sample());// * hamming[buffercount];//(outdata.channel[LEFT]);
		//X[CH1 + (t*N) + buffercount].numbers[REAL] = (float)outdata.channel[LEFT];// * hamming[buffercount];//(outdata.channel[LEFT]);
		//X[CH1 + (t*N) + buffercount].numbers[IMAG] = 0.0;
		//X[CH2 + (t*N) + buffercount].numbers[REAL] = (float)outdata.channel[RIGHT];// * hamming[buffercount];//(outdata.channel[LEFT]);
		//X[CH2 + (t*N) + buffercount++].numbers[IMAG] = 0.0;

		x[CH1_T + (t*N_INT) + buffercount] = (float)outdata.channel[RIGHT];
		x[CH2_T + (t*N_INT) + buffercount++] = (float)outdata.channel[LEFT];
	}
	
	if (buffercount >= N_INT)
	{
		buffercount = 0;
    	//bufferfull = 1;
		t++; // Increase time index
	}
}

     
void main(void)
{	
	unsigned int index=0, n=0, m=0;
	unsigned short k=0; // k is the freqency bin index
	complexpair *w_ptr = (complexpair*)w;	
	union complexdata *X1_ptr=&X[CH1], *X2_ptr=&X[CH2]; // Pointers to each individual channel

	
	// PCA variables
	COMPLEX mean1, mean2, temp[2]; 
	COMPLEX d[2], E[2][2], Rxx[2][2];
	COMPLEX_DBL Rxx_dbl[2][2];
	COMPLEX dbl_conver; //dbl_conver[2] is used to typecast type of COMPLEX to COMPLEX_DBL
	COMPLEX W_temp[4],W_inv[4],Q_inv[4];
	float D[2]={0.0, 0.0};
	
	// Test variables
/*	float test = 0.0;//, a = 2.1, b = 1.3;
	float test_r=0.0,test_i=0.0;
	COMPLEX c,e;
	e.real=0.4;
	e.imag=0.1;
	c.real=1.1;
	c.imag=2.1;
	cmplx_mult_add(c, e, e, e, &test_r, &test_i);
	temp[1] = cmplx_mult(c, e);
	*/
	//short brev[8];
	unsigned int prevGIE; 
	#pragma DATA_ALIGN(test, 8)
	float test[32] = {0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	#pragma DATA_ALIGN(w8, 8)
	double w8[16]={1.0000, 0, 0.7071,0.7071,0.0000,1.0000,-0.7071,0.7071,-1.0000, 0.0000,-0.7071,-0.7071,-0.0000,-1.0000,0.7071,-0.7071};
	float y[16];
	#pragma DATA_ALIGN(brev,8)
	/*short brev[64] = {
		0x0, 0x20, 0x10, 0x30, 0x8, 0x28, 0x18, 0x38,
		0x4, 0x24, 0x14, 0x34, 0xc, 0x2c, 0x1c, 0x3c,
		0x2, 0x22, 0x12, 0x32, 0xa, 0x2a, 0x1a, 0x3a,
		0x6, 0x26, 0x16, 0x36, 0xe, 0x2e, 0x1e, 0x3e,
		0x1, 0x21, 0x11, 0x31, 0x9, 0x29, 0x19, 0x39,
		0x5, 0x25, 0x15, 0x35, 0xd, 0x2d, 0x1d, 0x3d,
		0x3, 0x23, 0x13, 0x33, 0xb, 0x2b, 0x1b, 0x3b,
		0x7, 0x27, 0x17, 0x37, 0xf, 0x2f, 0x1f, 0x3f
	};*/
	//short brev1024[1024];
	//short brev[8] = {0,0,0,0,0,0,0,0};
	short brev[8] = {0,4,2,6,1,5,3,7};
	//short brev[8] = {0,2,1,3,0,0,0,0};
	//short brev[16] = {0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};
	//bitrev_index(&brev1024[0], 1024) ;
	//DSPF_sp_bitrev_cplx((double*)w, &brev1024[0], 1024);
	
	Xstart_ptr = &(X[0].cart);
	comm_intr(); 

	DSK6713_LED_off(0); // Make sure the LEDs are off
	DSK6713_LED_off(1);
	DSK6713_LED_off(2);
	DSK6713_LED_off(3);	
	
	/*
	CACHE_enableCaching(CACHE_CE00);
  	CACHE_setL2Mode(CACHE_64KCACHE);
	
	// Invalidate L1D, L1P, and L2 cache    
    CACHE_wbInvL1d((void *)0x0, 65536, CACHE_WAIT);
    CACHE_invAllL1p();
    CACHE_wbInvAllL2(CACHE_WAIT);    
	


	tw_genSPxSPfft(&w8[0], 8);
	
	bitrev_index(brev,8);
	prevGIE = IRQ_globalDisable(); // Turn off global interrupts
		
		// Bit reverse  frequency domain data
	//DSPF_sp_bitrev_cplx((double*)test, brev, 8);//N_INT=1024
 	DSPF_sp_ifftSPxSP(8, &test[0], &w8[0], y, brev, 2, 0, 8);
	IRQ_globalRestore(prevGIE);
  	*/
  	
  	
  	
  	
	for(n=STFT_SIZE;n>STFT_SIZE-N;n--) // write zeros on the end of the buffer 
	{									// Is this still necessary? 7-4-2013
		X[CH1 + n].cart.real=0.0;// 
		X[CH1 + n].cart.imag=0.0;// 
	}
	n=0;	



	
	while(t<TIME_BLOCKS_INT);// Input data loop
	
	

	m=0; // Don't think this a problem using m here?
	// 'Offline' STFT goes here
	for(n=0; n<((N_INT*TIME_BLOCKS_INT)-(N_INT/2)); n+=(N_INT/2)) // N/2 for 50% overlapping 
	{	
		
		// In order to implement the window you need to loop around every value and multiply it by the relevant coefficient 
		for(k=0;k<N_INT;k++)
		{

			buffer1[k].cart.real = hamming[k] * x[CH1_T + n+k];
			buffer2[k].cart.real = hamming[k] * x[CH2_T + n+k];
			buffer1[k].cart.imag = 0.0;
			buffer2[k].cart.imag = 0.0;
		}
		// BUT! memcpy seems to be more efficient - this way you can't use the window
		//memcpy(&buffer, &X[n], N*sizeof(complexpair)); // Copy full 1024 time domain points
		
		// Calculate 1024 point FFT on current buffers
		jack_fft(&buffer1->cart, N_INT, &w_ptr->cart);
		jack_fft(&buffer2->cart, N_INT, &w_ptr->cart);
		
		memcpy(&X[CH1 + n + m].cart, &buffer1[0].cart, N*sizeof(complexpair)); 
		memcpy(&X[CH2 + n + m].cart, &buffer2[0].cart, N*sizeof(complexpair));
		m++;
	}
	
	DSK6713_LED_on(1);
		
	
	
	/* PCA STARTS HERE - 2*2 case only*/
	for(k=0;k<N;k++)// Loop around half the number of frequency bins
	{
		D[0] = 0.0;
		D[1] = 0.0;
		d[0].real = 0.0;
		d[0].imag = 0.0;
		d[1].real = 0.0;
		d[1].imag = 0.0;
		E[0][0].real = 0.0;
		E[0][0].imag = 0.0;
		E[0][1].real = 0.0;
		E[0][1].imag = 0.0;
		E[1][0].real = 0.0;
		E[1][0].imag = 0.0;
		E[1][1].real = 0.0;
		E[1][1].imag = 0.0;
		
		mean1.real = 0.0;
		mean1.imag = 0.0;
		mean2.real = 0.0;
		mean2.imag = 0.0;
		
		// Find the mean value for both channels - in addition check for NaNs
		for(m=0; m<TIME_BLOCKS; m++)
		{
			index = N*m + k;
			
			//Comparisons involving NaN ALWAYS return FALSE, so invert by using != to check for NaNs
			if(X[CH1 + index].cart.real != X[CH1 + index].cart.real)
			{
				X[CH1 + index].cart.real = 0.0;
			}
			if(X[CH1 + index].cart.imag != X[CH1 + index].cart.imag)
			{
				X[CH1 + index].cart.imag = 0.0;
			}
			if(X[CH2 + index].cart.real != X[CH2 + index].cart.real)
			{
				X[CH2 + index].cart.real = 0.0;
			}
			if(X[CH2 + index].cart.imag != X[CH2 + index].cart.imag)
			{
				X[CH2 + index].cart.imag = 0.0;
			}
			
			mean1.real += ((*(X1_ptr+index)).cart.real);
			mean1.imag += ((*(X1_ptr+index)).cart.imag);
			mean2.real += ((*(X2_ptr+index)).cart.real);
			mean2.imag += ((*(X2_ptr+index)).cart.imag);
		
		//	mean1.real += *X1_ptr->cart->real;
		//	mean1.imag += (X1_ptr + N2*m + k)->cart->imag;
		//	mean2.real += (X2_ptr + N2*m + k)->cart->real;
		//	mean2.imag += (X2_ptr + N2*m + k)->cart->imag;
		
		}
		mean2.real = mean2.real/(float)TIME_BLOCKS; 
		mean2.imag = mean2.imag/(float)TIME_BLOCKS; 
		mean1.real = mean1.real/(float)TIME_BLOCKS; 
		mean1.imag = mean1.imag/(float)TIME_BLOCKS; 
	
		// Initialise covariance matrix for the current frequency bin
		Rxx[0][0].real = 0.0;
		Rxx[0][0].imag = 0.0;
		Rxx[0][1].real = 0.0;
		Rxx[0][1].imag = 0.0;
		Rxx[1][0].real = 0.0;
		Rxx[1][0].imag = 0.0;
		Rxx[1][1].real = 0.0;
		Rxx[1][1].imag = 0.0;
		Rxx_dbl[0][0].real = 0.0;
		Rxx_dbl[0][0].imag = 0.0;
		Rxx_dbl[0][1].real = 0.0;
		Rxx_dbl[0][1].imag = 0.0;
		Rxx_dbl[1][0].real = 0.0;
		Rxx_dbl[1][0].imag = 0.0;
		Rxx_dbl[1][1].real = 0.0;
		Rxx_dbl[1][1].imag = 0.0;
		
		
		
		// Calculate covariance matrix
		for(m=0;m<TIME_BLOCKS;m++)// 2 by many matrix multiplied by many by 2 matrix
		{
			index = N*m + k;
			temp[0].real = 0.0;
			temp[0].imag = 0.0;
			temp[1].real = 0.0;
			temp[1].imag = 0.0;
			temp[0].real = X[CH1 + index].cart.real - mean1.real; // N2 as half the FFT data is thrown away (and the num of time blocks has doubled)
			temp[0].imag = X[CH1 + index].cart.imag - mean1.imag;
			temp[1].real = X[CH2 + index].cart.real - mean2.real;
			temp[1].imag = X[CH2 + index].cart.imag - mean2.imag;
			
			dbl_conver = cmplx_mult(temp[0], conj(temp[0]));// Take the conjugate here as cov(X) = XX^H / N when complex
			Rxx_dbl[0][0].real += (double)dbl_conver.real;
			Rxx_dbl[0][0].imag += (double)dbl_conver.imag;
			
			
			dbl_conver = cmplx_mult(temp[0], conj(temp[1]));// Take the conjugate here as cov(X) = XX^H / N when complex
			Rxx_dbl[0][1].real += (double)dbl_conver.real;
			Rxx_dbl[0][1].imag += (double)dbl_conver.imag;
			
			dbl_conver = cmplx_mult(temp[1], conj(temp[1]));// Take the conjugate here as cov(X) = XX^H / N when complex
			Rxx_dbl[1][1].real += (double)dbl_conver.real;
			Rxx_dbl[1][1].imag += (double)dbl_conver.imag;
		}
			
		Rxx[0][0].real = (float)(Rxx_dbl[0][0].real / (float)TIME_BLOCKS);	
		Rxx[0][0].imag = (float)(Rxx_dbl[0][0].imag / (float)TIME_BLOCKS);	
		Rxx[0][1].real = (float)(Rxx_dbl[0][1].real / (float)TIME_BLOCKS);
		Rxx[0][1].imag = (float)(Rxx_dbl[0][1].imag / (float)TIME_BLOCKS);	
		Rxx[1][1].real = (float)(Rxx_dbl[1][1].real / (float)TIME_BLOCKS);
		Rxx[1][1].imag = (float)(Rxx_dbl[1][1].imag / (float)TIME_BLOCKS);	
	
	
	
		Rxx[1][0] = Rxx[0][1]; // Covariance matrix is symetrical
		
		eig((&Rxx[0][0]), (&d[0]), (&E[0][0]));
	
			
		if(d[1].real>d[0].real)// This has replaced the sort function
		{// Sort in decending order
			// Swap the real part of the eigenvalues
			temp[0].real=d[0].real;
			d[0].real = d[1].real;
			d[1].real = temp[0].real;
			// Swap the imaginary part of the eigenvalues
			temp[0].imag=d[0].imag;
			d[0].imag = d[1].imag;
			d[1].imag = temp[0].imag;
			// Swap the eigenvectors
			temp[0].real = E[0][1].real;
			temp[0].imag = E[0][1].imag;
			temp[1].real = E[1][1].real;
			temp[1].imag = E[1][1].imag;
			E[0][1].real = E[0][0].real; // Put first column into the second
			E[0][1].imag = E[0][0].imag;
			E[1][1].real = E[1][0].real;
			E[1][1].imag = E[1][0].imag;
			E[0][0].real = temp[0].real;
			E[0][0].imag = temp[0].imag;
			E[1][0].real = temp[1].real; // Put temporary column into the first column to complete the swap
			E[1][0].imag = temp[1].imag;
		}
		//else - No need for an else condition - leave as is.
	
	// --- OLD WAY OF FINDING D
		// Find D^(-.5), though it's complex so we use de Moivre's formula
		//r = mag(d[0]); // magnitude of d[0]
		//theta = arg(d[0]); // argument of d[0]
		//D[0] = -(1/sqrt(r))*cos(theta/2); // only calculate the real part
		//D[0].imag = sqrt(r)*sin(-theta/2); // Imaginary part is dropped in original code
		
		//r = mag(d[1]);//sqrt(pow(d[1],2) + pow(d[1],2)); //pythagoras
		//theta = arg(d[1]);//atan(d[1].imag/d[1].real);
		//D[1] = -(1/sqrt(r))*cos(theta/2); // only calculate the real part
		//D[1].imag = sqrt(r)*sin(-theta/2);

	// -- NEW WAY -- Throw away imaginary part completely
		//D[0] = pow(d[0].real,-0.5);
		//D[1] = pow(d[1].real,-0.5);
//		D[0] = 1.0/sqrt(d[0].real);
//		D[1] = 1.0/sqrt(d[1].real);
		D[0] = FastInvSqrt(d[0].real); // Google fast inverse square root for more info
		D[1] = FastInvSqrt(d[1].real);
		
	
	
		// 2*2 matrix complex multiplication where the non-diagonal elements are zero - note D is real
		Q[4*k + 0].real = D[0] * E[0][0].real; 
		Q[4*k + 0].imag = D[0] * E[0][0].imag * (-1);// Hermitian transpose
		Q[4*k + 1].real = D[0] * E[1][0].real; 
		Q[4*k + 1].imag = D[0] * E[1][0].imag * (-1);
		Q[4*k + 2].real = D[1] * E[0][1].real; 
		Q[4*k + 2].imag = D[1] * E[0][1].imag * (-1);
		Q[4*k + 3].real = D[1] * E[1][1].real; 
		Q[4*k + 3].imag = D[1] * E[1][1].imag * (-1);
	
			
			
		// Remove mean from X - mean values are worked out 
		for(m=0; m<TIME_BLOCKS; m++)// 2 by many matrix multiplied by many by 2 matrix
		{
			index = N*m + k;
			X[CH1 + index].cart.real -= mean1.real;
			X[CH1 + index].cart.imag -= mean1.imag;
			X[CH2 + index].cart.real -= mean2.real;
			X[CH2 + index].cart.imag -= mean2.imag;
		}	
		
		// Applies the whitening matrix
		// Xp(:,:,k) = Q(:,:,k)*(X(:,:,k)-Xmean); <- MATLAB code for the following loop		
		for(m=0; m<TIME_BLOCKS; m++)// 2 by many matrix multiplied by many by 2 matrix
		{
			index = N*m + k;
			temp[0] = cmplx_add(cmplx_mult(Q[4*k + 0],X[CH1 + index].cart), cmplx_mult(Q[4*k + 1],X[CH2 + index].cart));
			temp[1] = cmplx_add(cmplx_mult(Q[4*k + 2],X[CH1 + index].cart), cmplx_mult(Q[4*k + 3],X[CH2 + index].cart));
			
			X[CH1 + index].cart = temp[0]; 
			X[CH2 + index].cart = temp[1];
		}
	
		Wp[4*k + 0].real = 1;// Intialise unmixing matrix at each frequency bin 
		Wp[4*k + 0].imag = 0;
		Wp[4*k + 1].real = 0;
		Wp[4*k + 1].imag = 0;
		Wp[4*k + 2].real = 0;
		Wp[4*k + 2].imag = 0;
		Wp[4*k + 3].real = 1;
		Wp[4*k + 3].imag = 0;	
	}
	DSK6713_LED_on(2);
	istft(&Xstart_ptr[0], x_rec, N, 40960, N);// Resynthesise the first source
	
	iva(&Xstart_ptr[0], Wp, N);// IVA algorithm in a separate function

	// Now convert Wp to the actual unmixing matrix W (unwhitening stage)
	for(k=0;k<N;k++)
	{
		index=4*k;
		W_temp[0] = cmplx_add(cmplx_mult(Wp[index + 0], Q[index + 0]), cmplx_mult(Wp[index + 1], Q[index + 2]));// Intialise unmixing matrix at each frequency bin 
		W_temp[1] = cmplx_add(cmplx_mult(Wp[index + 0], Q[index + 1]), cmplx_mult(Wp[index + 1], Q[index + 3]));		
		W_temp[2] = cmplx_add(cmplx_mult(Wp[index + 2], Q[index + 0]), cmplx_mult(Wp[index + 3], Q[index + 2]));
		W_temp[3] = cmplx_add(cmplx_mult(Wp[index + 2], Q[index + 1]), cmplx_mult(Wp[index + 3], Q[index + 3]));		
		
		
		inv_2x2(&W_temp[0], &W_inv[0]);
		// 2*2 matrix complex multiplication where the non-diagonal elements are zero
		Wp[index + 0] = cmplx_mult(W_inv[0], W_temp[0]);
		Wp[index + 1] = cmplx_mult(W_inv[0], W_temp[1]); 
		Wp[index + 2] = cmplx_mult(W_inv[1], W_temp[2]); 
		Wp[index + 3] = cmplx_mult(W_inv[1], W_temp[3]); 
	}	
	DSK6713_LED_on(3);	
	
	
	// Now recover original sources in the frequency by multiplying by the inverse of the whitening matrix
	for(k=0;k<N;k++)
	{	
		inv_2x2(&Q[4*k], &Q_inv[0]); // Work out the inverse of Q (whitening matrix at each frequency bin)
		
		for(m=0; m<TIME_BLOCKS; m++)// 2 by many matrix multiplied by many by 2 matrix
		{
			index = N*m + k;
			cmplx_mult_add(Q_inv[0], X[CH1 + index].cart, Q_inv[1], X[CH2 + index].cart, &temp[0].real, &temp[0].imag);
			cmplx_mult_add(Q_inv[2], X[CH1 + index].cart, Q_inv[3], X[CH2 + index].cart, &temp[1].real, &temp[1].imag);
			
			X[CH1 + index].cart = temp[0]; 
			X[CH2 + index].cart = temp[1];
		}
	}	
	
	// Unmixing filter multiplied the estimated source @ each freq needs to go here!
	
	
	
	istft(&Xstart_ptr[0], x, N, 40960, N);// Resynthesise the first source
	
		
	while(1);
}                                 //end of main()
