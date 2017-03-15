/* inline.c - URDME inline propensity functions. */

/* S. Engblom 2017-02-23 */

#include <stdlib.h>
#include "inline.h"

/*----------------------------------------------------------------------*/
double inlineProp(const int *xx,const double *k,const int *i,
		  const int *prS, size_t nS, double vol, int sd)
/* Inline propensity. Given constants k[0..2], indices i[0..2] and
   volume vol, this function computes an elementary propensity at
   state xx. */
{
  double res = 0.0;

  /* Check if this propensity is turned off. */
  if (sd == 0) return 0.0;
  for (size_t j = 0; j < nS; j++) if (sd == prS[j]) return 0.0;

  /* Main formula. */
  if (k[0] != 0.0) {
    if (i[0] != i[1]) res += k[0]*xx[i[0]]*xx[i[1]]/vol;
    else res += k[0]*xx[i[0]]*(xx[i[0]]-1)/(2.0*vol);
  }
  if (k[1] != 0.0) res += k[1]*xx[i[2]];

  res += k[2]*vol;
  return res;
}
/*----------------------------------------------------------------------*/
