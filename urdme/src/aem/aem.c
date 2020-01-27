/* aem.c - URDME AEM solver. */

/* S. Engblom 2019-11-12 (Nreplicas syntax) */
/* S. Engblom 2017-02-17 (Major revision, URDME 1.3, Comsol 5) */
/* P. Bauer and S. Engblom 2012-05-10 */

#include <string.h>
#include <math.h>
#include <float.h>

#include "inline.h"
#include "propensities.h"
#include "binheap.h"
#include "report.h"
#include "aem.h"

#if !defined(MALLOC) || !defined(FREE)
  #error "Must define MALLOC and FREE."
#endif

/* Global variables for debug statistics */
int updates = 0;
int wakeups = 0;
int react2diff = 0;
int react2react = 0;
int diff2react = 0;
int diff2diff = 0;

void calcTimes(double* time,double* infTime,double tt,double old_rate, 
	       double new_rate,seeds seed);

/*----------------------------------------------------------------------*/
void aem(const PropensityFun *rfun,
	 const size_t *irD,const size_t *jcD,const double *prD,
	 const int *u0,
	 const size_t *irN,const size_t *jcN,const int *prN,
	 const size_t *irG,const size_t *jcG,
	 const double *tspan,const size_t tlen,const size_t Nreplicas,
	 int *U,
	 const double *vol,const double *ldata,const double *gdata,
	 const int *sd,
	 const size_t Ncells,
	 const size_t Mspecies,const size_t Mreactions,
	 const size_t dsize,
	 int report_level,const unsigned *seed_uint,
	 const double *K,const int *I,
	 const size_t *jcS,const int *prS,const size_t M1)

/* Specification of the inputs, see nsm.c */
{
  double tt;
  double *rrate;
  double *Ddiag,*isdrate2;
  double *reactTimes,*diffTimes,*reactInf,*diffInf;

  double oldrate;

  int *reactNode,*reactHeap,*diffNode,*diffHeap,*xx;
  int *jcE,*irE,*jvec;

  int re,cnt,nSub;
  int diffHeapSize,reactHeapSize;

  int errcode = 0;
  long total_diffusion = 0,total_reactions = 0;

  int subvol,spec,to_vol,to_spec;
  size_t i,j,it,k;
  const size_t Ndofs = Ncells*Mspecies;

  seeds *reactSeeds,*diffSeeds;

  ReportFun report = &URDMEreportFun;

  /* Main allocation. */

  /* current state vector */
  xx = (int *)MALLOC(Ndofs*sizeof(int));

  /* Reaction rate matrix (Mreactions X Ncells). In rrate we store all
     propensities for chemical rections. */
  rrate = (double *)MALLOC(Mreactions*Ncells*sizeof(double));

  /* binary heap for reaction events */
  reactHeapSize = Ncells*Mreactions;
  reactTimes = (double *)MALLOC(reactHeapSize*sizeof(double));
  reactNode = (int *)MALLOC(reactHeapSize*sizeof(int));
  reactHeap = (int *)MALLOC(reactHeapSize*sizeof(int));
  reactSeeds = (seeds *)MALLOC(reactHeapSize*sizeof(seeds));

  /* pre-infinity times */
  reactInf = (double *)MALLOC(reactHeapSize*sizeof(double));

  /* fix to run diffusion-only tests */
  if (reactHeapSize == 0) {
    reactTimes = (double *)MALLOC(sizeof(double));
    reactTimes[0] = INFINITY;
  }

  /* The diagonal value of the D-matrix is used frequently. For
     efficiency, we store the negative of D's diagonal in Ddiag. */
  Ddiag = (double *)MALLOC(Ndofs*sizeof(double));
  for (i = 0,nSub = 0; i < Ndofs; i++) {
    Ddiag[i] = 0.0;
    for (j = jcD[i]; j < jcD[i+1]; j++)
      if (irD[j] == i) {
	Ddiag[i] = -prD[j];
	break;
      }
    if (Ddiag[i] > 0) nSub++;
  }

  /* diffHeapSize for unidirectional diffusion problems -> subtract
     non-zero'd Nrows only */
  diffHeapSize = jcD[Ndofs]-nSub;

  /* allocate space for diffusion-event heap */
  if (diffHeapSize < 0) diffHeapSize = 0;
  diffTimes = (double *)MALLOC(diffHeapSize*sizeof(double));
  diffNode = (int *)MALLOC(diffHeapSize*sizeof(int));
  diffHeap = (int *)MALLOC(diffHeapSize*sizeof(int));
  diffSeeds = (seeds *)MALLOC(diffHeapSize*sizeof(seeds));
  isdrate2 = (double *)MALLOC(diffHeapSize*sizeof(double));

  /* pre-infinity times */
  diffInf = (double *)MALLOC(diffHeapSize*sizeof(double));

  /* sparse matrix to keep track of diffusion events */
  jcE = (int *)MALLOC((Ndofs+1)*sizeof(int));
  irE = (int *)MALLOC(diffHeapSize*sizeof(int));
  jvec = (int *)MALLOC(diffHeapSize*sizeof(int));

  /* fix to run reaction-only tests */
  if (diffHeapSize == 0) {
    diffTimes = (double *)MALLOC(sizeof(double));
    diffTimes[0] = INFINITY;
  }

  if (report_level == 2) {
    PRINTF("Reaction heap size: %d\n",(int)reactHeapSize);
    PRINTF("Diffusion heap size: %d\n",(int)diffHeapSize);
  }

  /* Loop over Nreplicas cases. */
  for (k = 0; k < Nreplicas; k++) {
    it = 0;
    tt = tspan[0];
    
    /* Set xx to the initial state. */
    memcpy(xx,&u0[k*Ndofs],Ndofs*sizeof(int));

    /* set new master seed */
    srand(seed_uint[k]);

    /* initialize erand48 seeds: reactions... */
    for (i = 0; i < reactHeapSize; i++)
      for (j = 0; j < 3; j++)
	reactSeeds[i][j] = rand()%65535;

    /* diffusions */
    for (i = 0; i < diffHeapSize; i++)
      for (j = 0; j < 3; j++)
	diffSeeds[i][j] = rand()%65535;
    
    /* Calculate the propensity for every reaction and every
       subvolume. */
    for (i = 0; i < Ncells; i++) {
      for (j = 0; j < M1; j++)
	rrate[i*Mreactions+j] = 
	  inlineProp(&xx[i*Mspecies],&K[j*3],&I[j*3],&prS[jcS[j]],
		     jcS[j+1]-jcS[j],vol[i],sd[i]);
      for (; j < Mreactions; j++)
	rrate[i*Mreactions+j] = (*rfun[j-M1])(&xx[i*Mspecies],tt,vol[i],
					      &ldata[i*dsize],gdata,sd[i]);
    }

    /* times to next reaction event */
    for (i = 0; i < reactHeapSize; i++) {
      reactTimes[i] = -log(1.0-erand48(reactSeeds[i]))/rrate[i]+tspan[0];

      if (reactTimes[i] <= 0.0)
	reactTimes[i] = INFINITY;

      reactHeap[i] = reactNode[i] = i;
    }

    /* pre-infinity times */
    for (i = 0; i < reactHeapSize; i++) reactInf[i] = 0;

    /* reaction heap */
    initialize_heap(reactTimes,reactNode,reactHeap,reactHeapSize);

    /* queue all non-diagonal entries of D as diffusion events */
    for (cnt = 0, i = 0; i < Ndofs; i++) {
      jcE[i] = cnt;
      for (j = jcD[i]; j < jcD[i + 1]; j++)
	if (irD[j] != i) {
	  diffNode[cnt] = diffHeap[cnt] = cnt;
	  isdrate2[cnt] = xx[i]*prD[j];
	  diffTimes[cnt] = -log(1.0-erand48(diffSeeds[cnt]))/isdrate2[cnt]
	    +tspan[0];

	  if (diffTimes[cnt] <= 0.0)
	    diffTimes[cnt] = INFINITY;
	  irE[cnt] = irD[j];
	  jvec[cnt] = j;

	  cnt++;
	}
    }
    jcE[i] = cnt;

    /* pre-infinity times */
    for (i = 0; i < diffHeapSize; i++) diffInf[i] = 0;

    /* diffusion heap */
    initialize_heap(diffTimes,diffNode,diffHeap,diffHeapSize);

    /* Main loop. */
    for ( ; ;) {

      /* determine next event */
      if (reactTimes[0] <= diffTimes[0])
	tt = reactTimes[0];
      else
	tt = diffTimes[0];

      /* report data */
      if (tt >= tspan[it] || isinf(tt)) {
	for (; it < tlen && (tt >= tspan[it] || isinf(tt)); it++) {
	  if (report_level)
	    /* virtual time in [0,Nreplicas*(tspan[tlen-1]-tspan[0])] */
	    report(tspan[it]-tspan[0]+k*(tspan[tlen-1]-tspan[0]),
		   0.0,Nreplicas*(tspan[tlen-1]-tspan[0]),
		   total_diffusion,total_reactions,errcode,
		   report_level);
	  memcpy(&U[k*Ndofs*tlen+Ndofs*it],xx,Ndofs*sizeof(int));
	}
	if (it >= tlen) break; /* main exit */
      }

      if (reactTimes[0] <= diffTimes[0]) {
	/* Reaction event. */

	/* a) extract subvolume and reaction re */
	subvol = reactNode[0]/Mreactions;
	re = reactNode[0]%Mreactions;

	/* b) Update the state of the subvolume subvol and update
	   dependent diffusion events. */
	/* loop over species that are affected by the reaction */
	for (i = jcN[re]; i < jcN[re+1]; i++) {

	  /* state update */
	  xx[subvol*Mspecies+irN[i]] += prN[i];
	  if (xx[subvol*Mspecies+irN[i]] < 0) errcode = 1;

	  /* update all dependent diffusion events */
	  for (cnt = jcE[subvol*Mspecies+irN[i]];
	       cnt < jcE[subvol*Mspecies+irN[i]+1]; cnt++) {
	    /* update the diffusion rate for the affected species &
	       subvolume */
	    oldrate = isdrate2[cnt];
	    isdrate2[cnt] += prD[jvec[cnt]] * prN[i];

	    /* updated times and reorder the heap */
	    calcTimes(&diffTimes[diffHeap[cnt]],&diffInf[cnt],tt,oldrate,
		      isdrate2[cnt],diffSeeds[cnt]);
	    update(diffHeap[cnt],diffTimes,diffNode,diffHeap,diffHeapSize);
	    react2diff++;
	  }
	}

	/* c) update dependent reaction events */
	for (i = jcG[Mspecies+re]; i < jcG[Mspecies+re+1]; i++) {
	  j = irG[i];
	  if (j != re) { /* see code underneath */
	    oldrate = rrate[subvol*Mreactions+j];
	    if (j < M1)
	      rrate[subvol*Mreactions+j] = 
		inlineProp(&xx[subvol*Mspecies],
			   &K[j*3],&I[j*3],&prS[jcS[j]],
			   jcS[j+1]-jcS[j],vol[subvol],sd[subvol]);
	    else
	      rrate[subvol*Mreactions+j] =
		(*rfun[j-M1])(&xx[subvol*Mspecies],tt,
			      vol[subvol],&ldata[subvol*dsize],gdata,sd[subvol]);

	    /* update times and reorder the heap */
	    calcTimes(&reactTimes[reactHeap[subvol*Mreactions+j]],
		      &reactInf[subvol*Mreactions+j],tt,oldrate,
		      rrate[subvol*Mreactions+j],
		      reactSeeds[subvol*Mreactions+j]);
	    update(reactHeap[subvol*Mreactions+j],reactTimes,reactNode,
		   reactHeap,reactHeapSize);
	    react2react++;
	  }
	}
	/* finish with j = re (the one that just happened), which need
	   not be in the dependency graph but must be updated
	   nevertheless */
	j = re;
	oldrate = rrate[subvol*Mreactions+j];
	if (j < M1)
	  rrate[subvol*Mreactions+j] = 
	    inlineProp(&xx[subvol*Mspecies],
		       &K[j*3],&I[j*3],&prS[jcS[j]],
		       jcS[j+1]-jcS[j],vol[subvol],sd[subvol]);
	else
	  rrate[subvol*Mreactions+j] =
	    (*rfun[j-M1])(&xx[subvol*Mspecies],tt,
			  vol[subvol],&ldata[subvol*dsize],gdata,sd[subvol]);

	calcTimes(&reactTimes[reactHeap[subvol*Mreactions+j]],
		  &reactInf[subvol*Mreactions+j],tt,oldrate,
		  rrate[subvol*Mreactions+j],
		  reactSeeds[subvol*Mreactions+j]);
	update(reactHeap[subvol*Mreactions+j],reactTimes,reactNode,
	       reactHeap,reactHeapSize);
	react2react++;

	total_reactions++; /* counter */
      }
      else {
	/* Diffusion event. */

	/* determine parameters of the event */
	for (i = 0; diffNode[0] >= jcE[i]; i++)
	  ;

	subvol = (i-1)/Mspecies;
	spec = (i-1)%Mspecies;
	to_vol = irE[diffNode[0]]/Mspecies;
	to_spec = irE[diffNode[0]]%Mspecies;

	/* Execute the diffusion event (check for negative elements). */
	xx[subvol*Mspecies+spec]--;
	if (xx[subvol*Mspecies+spec] < 0) errcode = 2;
	xx[to_vol*Mspecies+to_spec]++;

	/* Recalculate the reaction rates using dependency graph G. */
	for (i = jcG[spec]; i < jcG[spec+1]; i++) {
	  j = irG[i];
	  oldrate = rrate[subvol*Mreactions+j];
	  if (j < M1)
	    rrate[subvol*Mreactions+j] = 
	      inlineProp(&xx[subvol*Mspecies],
			 &K[j*3],&I[j*3],
			 &prS[jcS[j]],jcS[j+1]-jcS[j],
			 vol[subvol],sd[subvol]);
	  else
	    rrate[subvol*Mreactions+j] =
	      (*rfun[j-M1])(&xx[subvol*Mspecies],tt,
			    vol[subvol],&ldata[subvol*dsize],gdata,sd[subvol]);

	  /* Update reaction waiting time in outgoing subvolume */
	  calcTimes(&reactTimes[reactHeap[subvol*Mreactions+j]],
		    &reactInf[subvol*Mreactions+j],tt,oldrate,
		    rrate[subvol*Mreactions+j],
		    reactSeeds[subvol*Mreactions+j]);
	  update(reactHeap[subvol*Mreactions+j],reactTimes,reactNode,
		 reactHeap,reactHeapSize);

	  oldrate = rrate[to_vol*Mreactions+j];
	  if (j < M1)
	    rrate[to_vol*Mreactions+j] = 
	      inlineProp(&xx[to_vol*Mspecies],
			 &K[j*3],&I[j*3],
			 &prS[jcS[j]],jcS[j+1]-jcS[j],
			 vol[to_vol],sd[to_vol]);
	  else
	    rrate[to_vol*Mreactions+j] =
	      (*rfun[j-M1])(&xx[to_vol*Mspecies],tt,
			    vol[to_vol],&ldata[to_vol*dsize],gdata,sd[to_vol]);

	  /* Update reaction waiting time in incoming subvolume */
	  calcTimes(&reactTimes[reactHeap[to_vol*Mreactions+j]],
		    &reactInf[to_vol*Mreactions+j],tt,oldrate,
		    rrate[to_vol*Mreactions+j],
		    reactSeeds[to_vol*Mreactions+j]);
	  update(reactHeap[to_vol*Mreactions+j],reactTimes,reactNode,
		 reactHeap,reactHeapSize);
	  diff2react += 2;
	}

	/* Update all diffusion events of affected species in the
	   outgoing subvolume */
	for (cnt = jcE[subvol * Mspecies + spec];
	     cnt < jcE[subvol * Mspecies + spec + 1]; cnt++) {
	  oldrate = isdrate2[cnt];
	  isdrate2[cnt] -= prD[jvec[cnt]];

	  calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], tt, oldrate,
		    isdrate2[cnt], diffSeeds[cnt]);

	  update(diffHeap[cnt], diffTimes, diffNode, diffHeap, diffHeapSize);
	  diff2diff++;
	}

	/* Update all diffusion events of affected species in the
	   incoming subvolume */
	for (cnt = jcE[to_vol * Mspecies + spec];
	     cnt < jcE[to_vol * Mspecies + spec + 1]; cnt++) {
	  oldrate = isdrate2[cnt];

	  isdrate2[cnt] += prD[jvec[cnt]];

	  calcTimes(&diffTimes[diffHeap[cnt]], &diffInf[cnt], tt, oldrate,
		    isdrate2[cnt], diffSeeds[cnt]);

	  update(diffHeap[cnt], diffTimes, diffNode, diffHeap, diffHeapSize);
	  diff2diff++;
	}

	total_diffusion++; /* counter */
      }

      if (errcode) {
	/* Report the error and exit. */
	memcpy(&U[k*Ndofs*tlen+Ndofs*it],xx,Ndofs*sizeof(int));
	report(tt-tspan[0]+k*(tspan[tlen-1]-tspan[0]),
	       0.0,Nreplicas*(tspan[tlen-1]-tspan[0]),
	       total_diffusion,total_reactions,errcode,
	       report_level);
	break;
      }
    
    } /* Main simulation loop: for ( ; ; ) */
  } /* Loop over Nreplicas cases: for (k = 0; k < Nreplicas; k++) */

  if (report_level == 2) {
    PRINTF("Updates: %d Wake-ups: %d \n",updates,wakeups);
    PRINTF("React2React: %d React2diff: %d \n",react2react,react2diff);
    PRINTF("Diff2diff: %d Diff2react: %d \n",diff2diff,diff2react);
  }

  FREE(jvec);
  FREE(irE);
  FREE(jcE);

  FREE(diffInf);
  FREE(isdrate2);
  FREE(diffSeeds);
  FREE(diffHeap);
  FREE(diffNode);
  FREE(diffTimes);

  FREE(Ddiag);

  FREE(reactInf);
  FREE(reactSeeds);
  FREE(reactHeap);
  FREE(reactNode);
  FREE(reactTimes);

  FREE(rrate);
  FREE(xx);
}
/*-----------------------------------------------------------------------*/
void calcTimes(double* time,double* infTime,double tt,double old_rate,
	       double new_rate,seeds seed)
/* Calculate update of waiting times, including sleeping times. */
{
  double oldtime = time[0];

  if (isinf(oldtime)) {
    if (infTime[0] == 0) // Waking up first time
      time[0] = -log(1.0-erand48(seed))/new_rate+tt;
    else { // Waking up the 2nd..nth time
      if (new_rate > 0.0)
	time[0] = tt+(infTime[0]/new_rate);
    }
    wakeups++;
  }
  else {
    if (new_rate >= DBL_MIN) {
      if (time[0] == tt) // Regular update of current event
	time[0] = -log(1.0-erand48(seed))/new_rate+tt;
      else
	// Regular update of dependent events (rescaling)
	time[0] = ((old_rate/new_rate)*(time[0]-tt))+tt;
      updates++;
    }
    else { // Next event time set to infinity
      infTime[0] = (oldtime-tt)*old_rate;
      time[0] = INFINITY;
    }
  }
}
/*-----------------------------------------------------------------------*/
