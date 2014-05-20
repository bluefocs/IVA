#ifndef FOURIER_H_
#define FOURIER_H_

// Function prototype
void istft(COMPLEX *X, float *xtime, int nfreq, int time_len, int overlap);
void stft(COMPLEX *X, float *xtime, int nfreq, int time_len, int overlap);

#endif /*FOURIER_H_*/
