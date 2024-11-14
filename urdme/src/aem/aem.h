/* aem.h - Header file for AEM-solver. */

/* S. Engblom 2024-05-10 (data_time, ldata_time, gdata_time) */
/* S. Engblom 2019-11-12 (Nreplicas syntax) */
/* S. Engblom 2017-02-17 (Major revision, URDME 1.3, Comsol 5) */
/* P. Bauer and S. Engblom 2012-05-10 */

#ifndef AEM_H
#define AEM_H

typedef unsigned short int seeds[3];

void aem(const PropensityFun *rfun,
	 const size_t *irD,const size_t *jcD,const double *prD,
	 const int *u0,
	 const size_t *irN,const size_t *jcN,const int *prN,
	 const size_t *irG,const size_t *jcG,
	 const double *tspan,const size_t tlen,const size_t Nreplicas,
	 int *U,
	 const double *vol,const double *ldata,const double *gdata,
	 const double *data_time,const double *ldata_time,const double *gdata_time,
	 const int *sd,const size_t Ncells,
	 const size_t Mspecies,const size_t Mreactions,
	 const size_t ldsize,const size_t ldtsize,const size_t gdtsize,
	 const size_t dtlen,
	 int report_level,const unsigned *seed_uint,
	 const double *K,const int *I,
	 const size_t *jcS,const int *prS,const size_t M1
	 );

#endif /* AEM_H */
