/* ssa.c - URDME SSA solver. */

/* S. Engblom 2024-05-07 (data_time, ldata_time, gdata_time) */
/* S. Engblom 2019-11-15 (Nreplicas, multiple seeds syntax) */
/* S. Engblom 2017-02-22 */

#include <string.h>
#include <math.h>

#include "propensities.h"
#include "inline.h"
#include "report.h"
#include "ssa.h"

#if !defined(MALLOC) || !defined(FREE)
  #error "Must define MALLOC and FREE."
#endif

/*----------------------------------------------------------------------*/
void ssa(const PropensityFun *rfun, 
	 const int *u0,
	 const size_t *irN,const size_t *jcN,const int *prN,
	 const size_t *irG,const size_t *jcG,
	 const double *tspan,const size_t tlen,const size_t Nreplicas,
	 int *U,
	 const double *vol,const double *ldata,const double *gdata,
	 const double *data_time,const double *ldata_time,const double *gdata_time,
	 const int *sd,
	 const size_t Ncells,
	 const size_t Mspecies,const size_t Mreactions,
	 const size_t ldsize,const size_t ldtsize,const size_t gdtsize,
	 const size_t dtlen,
	 int report_level,const long *seed_long,
	 const double *K,const int *I,
	 const size_t *jcS,const int *prS,const size_t M1
	 )
/* Specification of the inputs, see nsm.c */
{
  /* reaction rates and sum of rates */
  double srrate,*rrate;

  /* state vector */
  int *xx;

  /* stats */
  long int total_reactions = 0;
  int errcode = 0;
  const int event = -1; /* only reactions are handled by SSA */
  const size_t Ndofs = Ncells*Mspecies;

  /* reporter */
  ReportFun report = &URDMEreportFun;

  /* allocate state and rate vectors */
  xx = (int *)MALLOC(Mspecies*sizeof(int));
  rrate = (double *)MALLOC(Mreactions*sizeof(double));

  /* Loop over Nreplicas cases. */
  for (int k = 0; k < Nreplicas; k++) {
    /* set new master seed */
    srand48(seed_long[k]);

    /* main loop over the (independent) cells */
    for (size_t subvol = 0; subvol < Ncells; subvol++) {
      /* Set (xx,tt) to the initial state. */
      size_t it = 0;
      double tt = tspan[0];
      memcpy(xx,&u0[Mspecies*subvol+k*Ndofs],Mspecies*sizeof(int));

      /* Pointer into data_time, ldata_time, gdata_time. */
      size_t timeix = 0;

      /* Calculate the propensity for every reaction. Store the sum of
	 the reaction intensities in srrate. */
      size_t j;
      srrate = 0.0;
      for (j = 0; j < M1; j++) {
	rrate[j] = inlineProp(xx,&K[j*3],&I[j*3],&prS[jcS[j]],
			      jcS[j+1]-jcS[j],vol[subvol],sd[subvol]);
	srrate += rrate[j];
      }
      for (; j < Mreactions; j++) {
	rrate[j] = (*rfun[j-M1])(xx,tt,vol[subvol],
				 &ldata[subvol*ldsize],gdata,
				 &ldata_time[(timeix*Ncells+subvol)*ldtsize],
				 &gdata_time[timeix*gdtsize],
				 sd[subvol]);
	srrate += rrate[j];
      }

      /* Main simulation loop. */
      for ( ; ; ) {

	/* time for next reaction */
	double dt = -log(1.0-drand48())/srrate;

	/* rather increase timeix first? */
	if (timeix+1 < dtlen && tt+dt >= data_time[timeix+1]) {
	  tt = data_time[++timeix];
	  dt = -1.0; /* signal nil event */

	  /* recalculate srrate and rrate using dependency graph */
	  for (int i = jcG[Mspecies+Mreactions]; i < jcG[Mspecies+Mreactions+1]; i++) {
	    const int j = irG[i];
	    const double old = rrate[j];
	    if (j < M1)
	      srrate += (rrate[j] = 
			 inlineProp(xx,&K[j*3],&I[j*3],&prS[jcS[j]],
				    jcS[j+1]-jcS[j],
				    vol[subvol],sd[subvol]))-old;
	    else
	      srrate += (rrate[j] = 
			 (*rfun[j-M1])(xx,tt,vol[subvol],
				       &ldata[subvol*ldsize],gdata,
				       &ldata_time[(timeix*Ncells+subvol)*ldtsize],
				       &gdata_time[timeix*gdtsize],
				       sd[subvol]))-old;
	  }
	}
	else
	  tt += dt;

	/* Store solution if the global time counter tt has passed the
	   next time in tspan. */
	if (tt >= tspan[it] || isinf(tt)) {
	  for (; it < tlen && (tt >= tspan[it] || isinf(tt)); it++)
	    memcpy(&U[k*Ndofs*tlen+Mspecies*subvol+Ndofs*it],xx,Mspecies*sizeof(int));

	  /* If the simulation has reached the final time, continue to
	     next subvolume. */
	  if (it >= tlen) break;
	}

	/* catch nil event: done reporting, so go for next event */
	if (dt == -1.0) goto next_event;

	/* a) Determine the reaction re that did occur. */
	const double rand = drand48()*srrate;
	double cum;
	int re;
	for (re = 0, cum = rrate[0]; re < Mreactions && rand > cum;
	     cum += rrate[++re]);

	/* elaborate floating point fix: */
	if (re >= Mreactions) re = Mreactions-1;
	if (rrate[re] == 0.0) {
	  /* go backwards and try to find first nonzero reaction rate */
	  for ( ; re > 0 && rrate[re] == 0.0; re--);

	  /* No nonzero rate found, but a reaction was sampled. This can
	     happen due to floating point errors in the iterated
	     recalculated rates. */
	  if (rrate[re] == 0.0) {
	    /* nil event: zero out and move on */
	    srrate = 0.0;
	    goto next_event;
	  }
	}

	/* b) Update the state of the subvolume. */
	for (int i = jcN[re]; i < jcN[re+1]; i++) {
	  xx[irN[i]] += prN[i];
	  if (xx[irN[i]] < 0) errcode = 1;
	}

	/* c) Recalculate srrate and rrate using dependency graph. */
	for (int i = jcG[Mspecies+re]; i < jcG[Mspecies+re+1]; i++) {
	  const int j = irG[i];
	  const double old = rrate[j];
	  if (j < M1)
	    srrate += (rrate[j] = 
		       inlineProp(xx,&K[j*3],&I[j*3],&prS[jcS[j]],
				  jcS[j+1]-jcS[j],
				  vol[subvol],sd[subvol]))-old;
	  else
	    srrate += (rrate[j] = 
		       (*rfun[j-M1])(xx,tt,vol[subvol],
				     &ldata[subvol*ldsize],gdata,
				     &ldata_time[(timeix*Ncells+subvol)*ldtsize],
				     &gdata_time[timeix*gdtsize],
				     sd[subvol]))-old;
	}

	total_reactions++; /* counter */

      next_event:
	/* Check for error codes. */
	if (errcode) {
	  /* Report the error that occurred and exit. */
	  memcpy(&U[k*Ndofs*tlen+Mspecies*subvol+Ndofs*it],xx,Mspecies*sizeof(int));
	  report(k*Ncells+subvol,0,Nreplicas*Ncells,
		 0,total_reactions,errcode,report_level);
	  break;
	}
      }
      if (report_level)
	report(k*Ncells+subvol,0,Nreplicas*Ncells,
	       0,total_reactions,0,report_level);
    }
  }

  FREE(rrate);
  FREE(xx);
}
/*----------------------------------------------------------------------*/
