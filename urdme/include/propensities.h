/* propensities.h - Declaration of URDME propensity functions. */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* S. Engblom 2017-02-16 (Major revision, URDME 1.3, Comsol 5) */
/* B. Drawert 2012-09-08 (revision) */
/* A. Hellander 2012-06-05 (revision) */  
/* J. Cullhed 2008-06-18 */

#ifndef PROPENSITIES_H
#define PROPENSITIES_H

#include "mex.h"

/* decide on the type of the state (integer/double) */
#ifdef UDS_
/* make by the deterministic solver? */
typedef double URDMEstate_t;
#else
/* the default: */
typedef int URDMEstate_t;
#endif

/* propensity function definition */
typedef double (*PropensityFun)(const URDMEstate_t *xstate,
				double time,double vol,
				const double *ldata,const double *gdata,
				int sd);

/* allocation and deallocation of propensity list */
PropensityFun *ALLOC_propensities(size_t Mreactions);
void FREE_propensities(PropensityFun* ptr);

#endif /* PROPENSITIES_H */
