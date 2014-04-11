#ifndef IVA_H_
#define IVA_H_
#include "definitions.h"

extern COMPLEX S[N2 * 2 * TIME_BLOCKS_50PC];
void iva(COMPLEX *Xp, COMPLEX *Wp, unsigned short nfreq);

//extern COMPLEX S[N2 * 2 * TIME_BLOCKS_50PC]; // Estimated sources variable
#endif /*IVA_H_*/
