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

#define DSK6713_AIC23_INPUT_MIC 0x0015
#define DSK6713_AIC23_INPUT_LINE 0x0011
#define N 1024
#define STFT_SIZE 80896
#define TIME_BLOCKS 79 // Actual no. of time blocks
#define TIME_BLOCKS_INT 40// Initial number of time blocks (no overlapping)
#define NSOURCES 2
#define LEFT 0 
#define RIGHT 1
#define REAL 0
#define IMAG 1
#define CH1 0
#define CH2 80896
#define CH1_T 0
#define CH2_T 40960
#endif /*DEFINITIONS_H_*/
