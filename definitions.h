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
#define N2 512
#define STFT_SIZE 40960//40000
#define TIME_BLOCKS 40//39
#define TIME_BLOCKS_50PC 79// The amount of time blocks the processing is actually done over
#define NSOURCES 2
#define LEFT 0 
#define RIGHT 1
#define REAL 0
#define IMAG 1
#define CH1 0
#define CH2 40960
#endif /*DEFINITIONS_H_*/
