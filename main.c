/*
 * main.c - written by Jack Harris Feburary-March-April 2014, Loughborough, UK.
 * 
 * 
 * 
 * LED 1) Reads in approx 5 seconds of audio data in progress.
 * LED 2) Indicates 1024 point FFT / STFT is complete with 50% overlapping with hanning 
 * 			window (see note).
 * LED 3) PCA complete.
 * LED 4) IVA complete, ISTFT in progress. Though feedthrough audio continues via the interrupt.
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
 *  - Inverse STFT still needs to be implemented. 
 * 				- Intial implementation completed 25-4-2014, 
 * 				- further testing 9-5-2014 (time domain signal still 'noisy')
 * 				- Buffer creation rewritten. 13-5-2014
 *  - (Possibly) expand frequency bins from 512 to 1024 (don't throw away half the 
 * 		frequency bins to speed up processing) - completed 14-4-2014
 * 	- Implement a complex multiply function in assembley. completed - 26-4-2014.
 * 	- Multiply the unmixing matrices by the estimated sources at each frequency bin. - Done 12-5-2014
 *  - Need to rewrite the exit clause for the IVA algorithm (checking the cost function).
 *  - Write FastInvSqrt function in linear assembly. 
 */
#include "dsk6713.h"
#include "dsk6713_led.h"
#include "dsk6713_dip.h"
#include "DSK6713_AIC23.h"	//codec-DSK interface support
//#include "fft.h"
#include "additional_math.h"
#include "iva.h"
#include "definitions.h"
#include "twiddles.h"
//#include "window.h"
#include "fourier.h"
#include "c67fastMath.h"

Uint32 fs=DSK6713_AIC23_FREQ_8KHZ;	//set sampling rate
DSK6713_AIC23_CodecHandle hCodec; // codec handle declaration
Uint16 inputsource=DSK6713_AIC23_INPUT_LINE; // select input source


union complexdata X[ NSOURCES * STFT_SIZE ]; // Used to store time data, freq domain data and whitened data
union complexdata X_org[NSOURCES*STFT_SIZE]; // Original X data for to be separated at the end
float x[NSOURCES * TIME_BLOCKS_INT * N_INT];
//float x_rec[TIME_BLOCKS_INT * N_INT];
unsigned short buffercount = 0;            //number of new input samples in buffer 
static unsigned short t = 0;
complexpair *X_ptr[TIME_BLOCKS];
COMPLEX *Xstart_ptr;

//
unsigned short outputstate=0;
unsigned int source_count=0; 
// Buffers for the FFT
complexpair buffer1[N_INT],buffer2[N_INT];	// Buffers for the FFTs of length 1024

// IVA variables - placed here so that the code works
COMPLEX Q[N*4];// 2*2 matrix at each freq bin
COMPLEX Wp[N*4];
#pragma DATA_SECTION(X,".EXT_RAM")
#pragma DATA_SECTION(X_org,".EXT_RAM")
#pragma DATA_SECTION(Q,".EXT_RAM")
#pragma DATA_SECTION(Wp,".EXT_RAM")
#pragma DATA_SECTION(buffer1,".EXT_RAM")
#pragma DATA_SECTION(buffer2,".EXT_RAM")
#pragma DATA_SECTION(x,".EXT_RAM")
//#pragma DATA_SECTION(x_rec,".EXT_RAM")
//#pragma DATA_SECTION(DSPF_sp_icfftr2_dif,".EXT_RAM")

// Function prototypes for any assembly routines
//extern COMPLEX cmplx_mult_sp(COMPLEX x, COMPLEX y, float *real, float *imag);


interrupt void c_int11(void)      //ISR
{
  	union {Uint32 uint; short channel[2];} outdata;
  	
	DSK6713_LED_on(0);		
	
	if (outputstate==0)/// This is the default
	{
		outdata.uint = input_sample(); // 
		output_sample(outdata.uint);
	}
	else if(outputstate==1)
	{
		output_left_sample(5*(short)(x[source_count++]));
		if(source_count>CH2)
		{
			source_count=0;
		}
	}
	else//outputstate==2
	{
		output_left_sample(5*(short)(x[CH2+source_count++]));
		if(source_count>CH2)
		{
			source_count=0;
		}
	}
	
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
	//complexpair *w_ptr = (complexpair*)w;	
	union complexdata *X1_ptr=&X[CH1], *X2_ptr=&X[CH2]; // Pointers to each individual channel
	COMPLEX *S1_ptr=&S[CH1], *S2_ptr=&S[CH2];
	
	// PCA variables
	COMPLEX mean1, mean2, temp[2]; 
	COMPLEX_DBL temp_dbl[2];
	COMPLEX d[2], E[2][2], Rxx[2][2];
	COMPLEX_DBL Rxx_dbl[2][2], E_dbl[2][2], d_dbl[2];//, Q_temp[4];
	COMPLEX_DBL dbl_conver; //dbl_conver[2] is used to typecast type of COMPLEX to COMPLEX_DBL
	COMPLEX W_temp[4],W_inv[4];//
	float D[2] = {0.0, 0.0};
	const float inv_timeblocks = 1/(float)TIME_BLOCKS;
	
	COMPLEX_DBL W_temp_dbl[4], X_temp[2];
	
	// Test variables
/*	float test_r=0.0,test_i=0.0;
	COMPLEX c,e;
	e.real=0.4;
	e.imag=0.1;
	c.real=1.1;
	c.imag=2.1;
	cmplx_mult_add(c, e, e, e, &test_r, &test_i);
	temp[1] = cmplx_mult(c, e);
*/
	
	
	Xstart_ptr = &(X[0].cart);
	comm_intr(); 

	DSK6713_LED_off(0); // Make sure the LEDs are off
	DSK6713_LED_off(1);
	DSK6713_LED_off(2);
	DSK6713_LED_off(3);	
	 	
  	
  	
	for(n=STFT_SIZE;n>STFT_SIZE-N;n--) // write zeros on the end of the buffer 
	{									// Is this still necessary? 7-4-2013
		X[CH1 + n].cart.real=0.0;// 
		X[CH1 + n].cart.imag=0.0;// 
	}
	n=0;	



	
	while(t<TIME_BLOCKS_INT);// Input data loop
	
	
/*
	
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
	*/
	stft(&X1_ptr->cart, x, N, 40960, 3*N_INT/4); // First 'microphone'
	stft(&X2_ptr->cart, &x[CH2_T], N, 40960, 3*N_INT/4); // Second 'microphone'
	
	memcpy(&X_org[0], &X[0], (NSOURCES*STFT_SIZE)*sizeof(complexpair)); // Save orginal STFT 
	DSK6713_LED_on(1);
	
	//istft(&X1_ptr->cart, &x[CH1], N, 40960, 3*N_INT/4);	// This is here to test the function
	
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
		mean2.real = mean2.real * inv_timeblocks; 
		mean2.imag = mean2.imag * inv_timeblocks; 
		mean1.real = mean1.real * inv_timeblocks; 
		mean1.imag = mean1.imag * inv_timeblocks; 
	
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
			temp_dbl[0].real = 0.0;
			temp_dbl[0].imag = 0.0;
			temp_dbl[1].real = 0.0;
			temp_dbl[1].imag = 0.0;
			temp_dbl[0].real = X[CH1 + index].cart.real - mean1.real; // N2 as half the FFT data is thrown away (and the num of time blocks has doubled)
			temp_dbl[0].imag = X[CH1 + index].cart.imag - mean1.imag;
			temp_dbl[1].real = X[CH2 + index].cart.real - mean2.real;
			temp_dbl[1].imag = X[CH2 + index].cart.imag - mean2.imag;
			
			dbl_conver = cmplx_mult_dbl(temp_dbl[0], conj_dbl(temp_dbl[0]));// Take the conjugate here as cov(X) = XX^H / N when complex
			Rxx_dbl[0][0].real += (double)dbl_conver.real;
			Rxx_dbl[0][0].imag += (double)dbl_conver.imag;
			
			
			dbl_conver = cmplx_mult_dbl(temp_dbl[0], conj_dbl(temp_dbl[1]));// Take the conjugate here as cov(X) = XX^H / N when complex
			Rxx_dbl[0][1].real += (double)dbl_conver.real;
			Rxx_dbl[0][1].imag += (double)dbl_conver.imag;
			
			dbl_conver = cmplx_mult_dbl(temp_dbl[1], conj_dbl(temp_dbl[1]));// Take the conjugate here as cov(X) = XX^H / N when complex
			Rxx_dbl[1][1].real += (double)dbl_conver.real;
			Rxx_dbl[1][1].imag += (double)dbl_conver.imag;
		}
			
		Rxx_dbl[0][0].real = Rxx_dbl[0][0].real * (double)inv_timeblocks;	
		Rxx_dbl[0][0].imag = Rxx_dbl[0][0].imag * (double)inv_timeblocks;	
		Rxx_dbl[0][1].real = Rxx_dbl[0][1].real * (double)inv_timeblocks;
		Rxx_dbl[0][1].imag = Rxx_dbl[0][1].imag * (double)inv_timeblocks;	
		Rxx_dbl[1][1].real = Rxx_dbl[1][1].real * (double)inv_timeblocks;
		Rxx_dbl[1][1].imag = Rxx_dbl[1][1].imag * (double)inv_timeblocks;	
	
	
	
		//Rxx[1][0] = Rxx[0][1]; // Covariance matrix is symetrical
		//Rxx[1][0].imag = Rxx[1][0].imag*(-1); 
		
		
		Rxx_dbl[1][0] = Rxx_dbl[0][1]; // Covariance matrix is symetrical
		Rxx_dbl[1][0].imag = Rxx_dbl[1][0].imag*(-1); 
		
		eig_dbl(&Rxx_dbl[0][0], d_dbl, &E_dbl[0][0]);
		
		d[0].real = (float)d_dbl[0].real;
		d[0].imag = (float)d_dbl[0].imag;
		d[1].real = (float)d_dbl[1].real;
		d[1].imag = (float)d_dbl[1].imag;
		E[0][0].real = (float)E_dbl[0][0].real;
		E[0][0].imag = (float)E_dbl[0][0].imag;
		E[0][1].real = (float)E_dbl[0][1].real;
		E[0][1].imag = (float)E_dbl[0][1].imag;
		E[1][0].real = (float)E_dbl[1][0].real;
		E[1][0].imag = (float)E_dbl[1][0].imag;
		E[1][1].real = (float)E_dbl[1][1].real;
		E[1][1].imag = (float)E_dbl[1][1].imag;

		
		// Fix scaling of the eigenvectors at first freq bin (always seems to be wrong because of eigenvectors)
	/*	if((k==0) || (k==(N-1)))
		{
			E[0][0].real = E[0][0].real * -1;
			E[0][0].imag = E[0][0].imag * -1;
			E[0][1].real = E[0][1].real * -1;
			E[0][1].imag = E[0][1].imag * -1;
			E[1][0].imag = E[1][0].imag * -1;
			E[1][0].real = E[1][0].real * -1;
			E[1][1].imag = E[1][1].imag * -1;
			E[1][1].real = E[1][1].real * -1;
			E[0][0].imag = E[0][0].imag * -1;
		}*/
			
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
			temp[0].real = E[0][1].real; // 'save' the second column
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
//		D[0] = FastInvSqrt(d[0].real); // Google fast inverse square root for more info
//		D[1] = FastInvSqrt(d[1].real);
		D[0] = rsqrtf(d[0].real);
		D[1] = rsqrtf(d[1].real);
				
		Q[4*k + 0].real = D[0] * (float)E[0][0].real; 
		Q[4*k + 0].imag = D[0] * (float)E[0][0].imag * (-1);// Hermitian transpose
		Q[4*k + 1].real = D[0] * (float)E[1][0].real; 
		Q[4*k + 1].imag = D[0] * (float)E[1][0].imag * (-1);
		Q[4*k + 2].real = D[1] * (float)E[0][1].real; 
		Q[4*k + 2].imag = D[1] * (float)E[0][1].imag * (-1);
		Q[4*k + 3].real = D[1] * (float)E[1][1].real; 
		Q[4*k + 3].imag = D[1] * (float)E[1][1].imag * (-1);
			
			
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
	/******* PCA ENDS HERE *******/
	
	DSK6713_LED_on(2);
	//istft(&Xstart_ptr[0], x_rec, N, 40960, N);// Resynthesise the first mixture to test istft
	
	iva(&Xstart_ptr[0], Wp, N);// IVA algorithm in a separate function

	// Now convert Wp to the actual unmixing matrix W (correct scaling of unmixing filter coefficients)
	for(k=0;k<N;k++)
	{
		index=4*k;
		W_temp[0] = cmplx_add(cmplx_mult(Wp[index + 0], Q[index + 0]), cmplx_mult(Wp[index + 1], Q[index + 2]));// Intialise unmixing matrix at each frequency bin 
		W_temp[1] = cmplx_add(cmplx_mult(Wp[index + 0], Q[index + 1]), cmplx_mult(Wp[index + 1], Q[index + 3]));		
		W_temp[2] = cmplx_add(cmplx_mult(Wp[index + 2], Q[index + 0]), cmplx_mult(Wp[index + 3], Q[index + 2]));
		W_temp[3] = cmplx_add(cmplx_mult(Wp[index + 2], Q[index + 1]), cmplx_mult(Wp[index + 3], Q[index + 3]));		
		
		
		inv_2x2(&W_temp[0], &W_inv[0]);
		
		// 2*2 matrix complex multiplication where the non-diagonal elements are zero - (discard W_inv[1] and W_inv[2])
		Wp[index + 0] = cmplx_mult(W_inv[0], W_temp[0]);
		Wp[index + 1] = cmplx_mult(W_inv[0], W_temp[1]); 
		Wp[index + 2] = cmplx_mult(W_inv[3], W_temp[2]); // Used to be W_inv[1], not sure why?
		Wp[index + 3] = cmplx_mult(W_inv[3], W_temp[3]); // Used to be W_inv[1], not sure why?
	}	
	DSK6713_LED_on(3);	
	
		
	// Now recover original sources in the frequency by multiplying by the scaled unmixing matrix
	for(k=0;k<N;k++)
	{		
		for(m=0; m<TIME_BLOCKS; m++)// 2 by many matrix multiplied by many by 2 matrix (in fact the multication is W*X at each freq bin)
		{
			index = N*m + k;
			//cmplx_mult_add(Wp[4*k + 0], X_org[CH1 + index].cart, Wp[4*k + 1], X_org[CH2 + index].cart, &temp[0].real, &temp[0].imag);
			//cmplx_mult_add(Wp[4*k + 2], X_org[CH1 + index].cart, Wp[4*k + 3], X_org[CH2 + index].cart, &temp[1].real, &temp[1].imag);
			
			
			W_temp_dbl[0].real = (double)(Wp[4*k + 0].real);
			W_temp_dbl[0].imag = (double)Wp[4*k + 0].imag;
			W_temp_dbl[1].real = (double)Wp[4*k + 1].real;
			W_temp_dbl[1].imag = (double)Wp[4*k + 1].imag;
			W_temp_dbl[2].real = (double)Wp[4*k + 2].real;
			W_temp_dbl[2].imag = (double)Wp[4*k + 2].imag;
			W_temp_dbl[3].real = (double)Wp[4*k + 3].real;
			W_temp_dbl[3].imag = (double)Wp[4*k + 3].imag;
			X_temp[0].real = (double)X_org[CH1 + index].cart.real;
			X_temp[0].imag = (double)X_org[CH1 + index].cart.imag;
			X_temp[1].real = (double)X_org[CH2 + index].cart.real;
			X_temp[1].imag = (double)X_org[CH2 + index].cart.imag;
			
			/* Seems that the filtering effect is inftroduced at this point due to the four complex multiplications above */
			temp_dbl[0] =cmplx_mult_dbl(W_temp_dbl[0], X_temp[0]);
			temp_dbl[1] =cmplx_mult_dbl(W_temp_dbl[1], X_temp[1]);
			S[CH1 + index].real = temp_dbl[0].real + temp_dbl[1].real;
			S[CH1 + index].imag = temp_dbl[0].imag + temp_dbl[1].imag;
			 
			temp_dbl[0] =cmplx_mult_dbl(W_temp_dbl[2], X_temp[0]);
			temp_dbl[1] =cmplx_mult_dbl(W_temp_dbl[3], X_temp[1]);			
			S[CH2 + index].real = temp_dbl[0].real + temp_dbl[1].real;
			S[CH2 + index].imag = temp_dbl[0].imag + temp_dbl[1].imag;
		}
	}	
	
	
	istft(S2_ptr, &x[CH2], N, 40960, 3*N_INT/4);// Resynthesise the second source
	istft(S1_ptr, &x[CH1], N, 40960, 3*N_INT/4);// Resynthesise the first source
	
		
	while(1)
	{
		if(DSK6713_DIP_get(1)==0)
		{
			outputstate=1;
		}
		else if(DSK6713_DIP_get(0)==0)	
		{
			outputstate=2;
		}
		else
		{
			outputstate=0;// Default - feedthrough
		}
	}
}                                 //end of main()
