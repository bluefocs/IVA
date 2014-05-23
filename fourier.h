#ifndef FOURIER_H_
#define FOURIER_H_

// Function prototype
void istft(COMPLEX *X, float *xtime, int nfreq, int time_len, int overlap);
void stft(COMPLEX *X, float *xtime, int nfreq, int time_len, int overlap);
void DSPF_sp_cfftr2_dit_c(float* x, float* w, short n); // C version of TI FFT
void DSPF_sp_icfftr2_dif_c(float* x, float * w, short n);
#endif /*FOURIER_H_*/
