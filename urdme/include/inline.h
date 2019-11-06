/* inline.h - URDME inline propensities. */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* S. Engblom 2017-02-23 */

#ifndef INLINE_H
#define INLINE_H

#include <stdlib.h>

#ifndef PROPENSITIES_H
#include "propensities.h"
#endif

double inlineProp(const URDMEstate_t *xx,const double *k,const int *i,
		  const int *prS, size_t nS, double vol, int sd);

#endif /* INLINE_H */
