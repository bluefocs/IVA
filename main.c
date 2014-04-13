/*
 * main.c - written by Jack Harris Feburary-March-April 2014, Loughborough, UK.
 * 
 * 
 * Interrupt - continuous 'feedthrough' audio system continuously taking audio data
 * from the line in and outputting it to the headphones.
 * 
 * LED 1) Reads in approx 5 seconds of audio data.
 * LED 2) Indicates 1024 point FFT / STFT is in progress with 50% overlapping and 
 * 			optional windowing (see note).
 * LED 3) PCA in progress.
 * LED 4) PCA is complete. Though feedthrough audio continues via the interrupt. IVA 
 * 			in progress 
 * 
 * Twiddle values for the FFT are generated in MATLAB (generate_twiddles.m) due to 
 * previous problems with the cos (cosine) function with the TI toolchain.
 * 
 * The STFT works by first assuming no overlap, the overlap is created by calculating 
 * each FFT in turn, but only writting the first half of the transformed signal back to 
 * memory (as the FFT of a real signal is symmetrical), then advancing half the FFT size
 * (in this case 1024/2), and repe
 * ating process to obtain the FFTs 'inbetween' the 
 * nonoverlapped FFTs.
 * 
 * The PCA and whitening has been fully debugged. Returns identity matrix when two 
 * corresponding frequency bins are compared (i.e. they are uncorrelated).
 * 
 * IVA has been been debugged. This has been verified with the output of the IVA 
 * algorithm in MATLAB (provided by Taesu Kim). The exit clause needs to be investigated,
 * this version runs too long/too short. However returns the correct Wp (or near enough).
 * 
 * The program can suffer from NaN propogating through each stage, this is POSSIBLY 
 * caused by the FFT code from Chassaing's book and dealing with large input values, 
 * in future versions scaling needs to be considered. Possibly turn down the volume 
 * of whatever is driving the audio in. 
 * 
 * Losely based on the original file in Rulph Chassaing's TMS320C6713 DSK book:
 * frames.c - basic illustration of triple-buffered
 *  N sample frame-based processing
 * 
 * TO DO LIST:
 *  - Inverse STFT still needs to be implemented.
 *  - (Possibly) expand frequency bins from 512 to 1024 (don't throw away half the 
 * 		frequency bins to speed up processing)
 * 
 */
#include "dsk6713.h"
#include "DSK6713_AIC23.h"	//codec-DSK interface support
#include "fft.h"
#include "additional_math.h"
#include "iva.h"
#include "definitions.h"
#include "twiddles.h"
#include "hamming.h"

Uint32 fs=DSK6713_AIC23_FREQ_8KHZ;	//set sampling rate
DSK6713_AIC23_CodecHandle hCodec; // codec handle declaration
Uint16 inputsource=DSK6713_AIC23_INPUT_LINE; // select input source


union complexdata X[NSOURCES*STFT_SIZE];// Used to store time data, freq domain data and whitened data
unsigned short buffercount = 0;            //number of new input samples in buffer 
unsigned short bufferfull = 0;                   //set by ISR to indicate iobuffer full
static unsigned short t = 0;
complexpair *X_ptr[TIME_BLOCKS];
COMPLEX *Xstart_ptr;
//complexpair** currentBlock_ptr = X_ptr; 


// IVA variables - placed here so that the code works
COMPLEX Q[N*4];// 2*2 matrix at each freq bin
COMPLEX Wp[N*4];
#pragma DATA_SECTION(X,".EXT_RAM")
#pragma DATA_SECTION(Q,".EXT_RAM")
#pragma DATA_SECTION(Wp,".EXT_RAM")

interrupt void c_int11(void)      //ISR
{
  	union {Uint32 uint; Uint16 channel[2];} outdata;
  	
	DSK6713_LED_on(0);		
	
	outdata.uint = input_sample(); // 
	output_sample(outdata.uint);
	
	//Original way of reading in data
	//(input_ptr + buffercount)->numbers[IMAG] = 0.0; 
	//(input_ptr + buffercount++)->numbers[REAL] = (float)(outdata.channel[LEFT]); 
		
	if(t<(TIME_BLOCKS))
	{	
		//(X_ptr[t] + buffercount)->numbers[IMAG] = 0.0; 
		//(X_ptr[t] + buffercount++)->numbers[REAL] = (float)(outdata.channel[LEFT]);	
		
		//X[CH1 + (t*N) + buffercount].numbers[REAL] = (float)(input_left_sample());// * hamming[buffercount];//(outdata.channel[LEFT]);
		X[CH1 + (t*N) + buffercount].numbers[REAL] = (float)outdata.channel[LEFT];// * hamming[buffercount];//(outdata.channel[LEFT]);
		X[CH1 + (t*N) + buffercount].numbers[IMAG] = 0.0;
		X[CH2 + (t*N) + buffercount].numbers[REAL] = (float)outdata.channel[RIGHT];// * hamming[buffercount];//(outdata.channel[LEFT]);
		X[CH2 + (t*N) + buffercount++].numbers[IMAG] = 0.0;
		
	}
	
	if (buffercount >= N)
	{
		buffercount = 0;
    	bufferfull = 1;
		t++; // Increase time index
	}
}

     
void main(void)
{	
	unsigned int index=0;
	unsigned short n=0,m=0,k=0; // k is the freqency bin index
	complexpair *w_ptr = (complexpair*)w;
	complexpair buffer1[N],buffer2[N];	// Buffers for the FFTs of length 1024
	union complexdata *X1_ptr=&X[CH1], *X2_ptr=&X[CH2]; // Pointers to each individual channel

	
	// PCA variables
	COMPLEX mean1, mean2, temp[2]; 
	COMPLEX d[2], E[2][2], Rxx[2][2];
	COMPLEX_DBL Rxx_dbl[2][2];
	COMPLEX dbl_conver; //dbl_conver[2] is used to typecast type of COMPLEX to COMPLEX_DBL
	COMPLEX W_temp[4],W_inv[4];
	float D[2]={0.0, 0.0};

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
	
	for(n=STFT_SIZE;n>STFT_SIZE-N;n--) // write zeros on the end of the buffer 
	{									// Is this still necessary? 7-4-2013
		X[CH1 + n].cart.real=0.0;// 
		X[CH1 + n].cart.imag=0.0;// 
	}
	n=0;	
	

	Xstart_ptr = &(X[0].cart);
	comm_intr(); 

	DSK6713_LED_off(0); // Make sure the LEDs are off
	DSK6713_LED_off(1);
	DSK6713_LED_off(2);
	DSK6713_LED_off(3);	
	
	
	while(t<TIME_BLOCKS)//was(t<TIME_BLOCKS-1)                        //frame processing loop
	{
		//while(bufferfull==0); //wait until buffer is full 
    	//bufferfull = 0;
    	
		/* FFTs would have gone here - now done offline 
		 * if any overlapping is required it'll need to go here. */
	}                               //end of first while
	DSK6713_LED_on(1);
	
	
	
	// 'Offline' STFT goes here
	//for(n=0; n<N*(TIME_BLOCKS-1); n+=N2)// N2 is half of N
	//for(n=0; n<(N2*TIME_BLOCKS_50PC); n+=N2)// N2 is half of N
	for(n=0; n<(N*TIME_BLOCKS); n+=N) // No overlapping 
	{	
		// In order to implement the window you need to loop around every value and multiply it by the relevant coefficient 
		for(k=0;k<N;k++)
		{
			buffer1[k].cart.real = hamming[k] * X[CH1 + n+k].cart.real;
			buffer2[k].cart.real = hamming[k] * X[CH2 + n+k].cart.real;
			buffer1[k].cart.imag = 0.0;
			buffer2[k].cart.imag = 0.0;
		}
		// BUT! memcpy seems to be more efficient - this way you can't use the window
		//memcpy(&buffer, &X[n], N*sizeof(complexpair)); // Copy full 1024 time domain points
		
		// Perform FFT on current buffers
		jack_fft(&buffer1->cart, N, &w_ptr->cart);// for CH1
		jack_fft(&buffer2->cart, N, &w_ptr->cart);// for CH2
		
		// Then copy 512 freq domain points back - i.e. store the first half of the FFT
		//memcpy(&X[CH1 + n], &buffer1, N2*sizeof(complexpair)); 
		//memcpy(&X[CH2 + n], &buffer2, N2*sizeof(complexpair)); 
		memcpy(&X[CH1 + n], &buffer1, N*sizeof(complexpair)); 
		memcpy(&X[CH2 + n], &buffer2, N*sizeof(complexpair)); 
	}
	
	/*// Old non-overlapping way of doing the STFT
	for(n=0;n<(TIME_BLOCKS-1);n++)
	{
		//fft(&X_ptr[n]->cart, N, &w_ptr->cart); // use 'n' as the time block index, so it doesn't clash with t
		jack_fft(&X_ptr[n]->cart, N, &w_ptr->cart);
	}*/
	DSK6713_LED_on(2);	
	
	
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
	
	DSK6713_LED_on(3);	
	
	iva(&Xstart_ptr[0], Wp, N);// IVA algorithm in a separate function

	// Now convert Wp to the actual unmixing matrix W (unwhitening stage?)
	for(k=0;k<N;k++)
	{
		W_temp[0] = cmplx_add(cmplx_mult(Wp[4*k + 0], Q[4*k + 0]), cmplx_mult(Wp[4*k + 1], Q[4*k + 2]));// Intialise unmixing matrix at each frequency bin 
		W_temp[1] = cmplx_add(cmplx_mult(Wp[4*k + 0], Q[4*k + 1]), cmplx_mult(Wp[4*k + 1], Q[4*k + 3]));		
		W_temp[2] = cmplx_add(cmplx_mult(Wp[4*k + 2], Q[4*k + 0]), cmplx_mult(Wp[4*k + 3], Q[4*k + 2]));
		W_temp[3] = cmplx_add(cmplx_mult(Wp[4*k + 2], Q[4*k + 1]), cmplx_mult(Wp[4*k + 3], Q[4*k + 3]));		
		
		
		inv_2x2(&W_temp[0], &W_inv[0]);
		// 2*2 matrix complex multiplication where the non-diagonal elements are zero
		Wp[4*k + 0] = cmplx_mult(W_inv[0], W_temp[0]);
		Wp[4*k + 1] = cmplx_mult(W_inv[0], W_temp[1]); 
		Wp[4*k + 2] = cmplx_mult(W_inv[1], W_temp[2]); 
		Wp[4*k + 3] = cmplx_mult(W_inv[1], W_temp[3]); 
	}			
	while(1);
}                                 //end of main()
