#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_
struct cmpx_dbl                       //complex data structure used by FFT
    {
    	double real;
    	double imag;
    };
typedef struct cmpx_dbl COMPLEX_DBL;

struct cmpx                       //complex data structure used by FFT
    {
    float real;
    float imag;
    };
typedef struct cmpx COMPLEX;

typedef union complexdata 
{
		float numbers[2];
		float number;
		long unsigned int full;
		COMPLEX cart;//Stands for cartesian - defined in fft.h
		
}complexpair;


// These definitions have changed for an STFT with a 3/4 overlap
#define DSK6713_AIC23_INPUT_MIC 0x0015
#define DSK6713_AIC23_INPUT_LINE 0x0011
#define N 513	// Number of frequency bins that the processing will be done across (N_INT/2) +1
#define N_INT 1024 // Original FFT length, INT = initial
#define STFT_SIZE 83106 // 513*162 // 513*156 //40527 // 513*79
#define TIME_BLOCKS 162//156 // Actual no. of time blocks
#define TIME_BLOCKS_INT 40// Initial number of time blocks (no overlapping)
#define NSOURCES 2
#define LEFT 0 
#define RIGHT 1
#define REAL 0
#define IMAG 1
#define CH1 0
#define CH2 83106//40527// This should be the same as STFT size
#define CH1_T 0
#define CH2_T 40960
#define PI 3.14159265358979323846
#endif /*DEFINITIONS_H_*/
