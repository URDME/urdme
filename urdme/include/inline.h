/* inline.h - URDME inline propensities. */

/* S. Engblom 2017-02-23 */

#ifndef INLINE_H
#define INLINE_H

#include <stdlib.h>

double inlineProp(const int *xx,const double *k,const int *i,
		  const int *prS, size_t nS, double vol, int sd);

#endif /* INLINE_H */
