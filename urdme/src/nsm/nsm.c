/* nsm.c - URDME NSM solver. */

/* S. Engblom 2024-05-08 (data_time, ldata_time, gdata_time) */
/* S. Engblom 2024-04-11 (generalized diffusion) */
/* S. Engblom 2019-11-12 (multiple seeds) */
/* S. Engblom 2018-02-10 (Nreplicas syntax) */
/* S. Engblom 2017-02-16 (Major revision, URDME 1.3, Comsol 5) */
/* S. Engblom 2014-06-10 (Revision) */
/* A. Hellander 2012-06-15 (Revision) */
/* P. Bauer and S. Engblom 2012-04-10 (Revision) */
/* A. Hellander 2009-11-24 (Revision) */
/* J. Cullhed 2008-06-18 (rdme.c) */

#include <string.h>
#include <math.h>

#include "propensities.h"
#include "binheap.h"
#include "inline.h"
#include "report.h"
#include "nsm.h"

#if !defined(MALLOC) || !defined(FREE)
  #error "Must define MALLOC and FREE."
#endif

/*----------------------------------------------------------------------*/
void nsm(const PropensityFun *rfun,
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
	 int report_level,const long *seed_long,
	 const double *K,const int *I,
	 const size_t *jcS,const int *prS,const size_t M1
	 )
/* Specification of the inputs:

Ncells
  Number of subvolumes.

Mspecies
  Number of species.

Hence Ndofs = Ncells*Mspecies.

Mreactions
  Total number of reactions.

M1.
  Number of inline propensities.

dsize
  Size of data vector sent to propensities.

tlen
  Number of sampling points in time.

Nreplicas
  Number of replicas, the 3rd dimension of u0.

report_level
  The desired degree of feedback during simulations. 0, 1, and 2 are
  currently supported options.

seed_long
  Vector of Nreplicas random seed values, passed to srand48.

Diffusion matrix D. Double sparse (Ndofs X Ndofs). 
  Macroscopic diffusion matrix. For i ~= j, D(i,j) is the linear
  transport rate from dof #j to dof #i. -D(j,j) is the total outward
  rate from dof #j.

Initial state vector u0. Integer (Mspecies X Ncells).
  Gives the initial copy number of the species in each subvolume.

Stochiometric matrix N. Integer sparse (Mspecies X Nreactions).
  N(:,j) describes how reaction j changes the number of species.

Dependency graph G. Integer sparse (Mreactions X Mspecies+Mreactions).
  G(i,Mspecies+j) is non-zero if executing reaction j means that
  reaction i needs to be re-evaluated. The first Mspecies columns of G
  similarily cover diffusion events.

Inline propensity coefficients K. Double (3 X M1)
  K(:,i) gives the coefficients to reaction #i. Where i is in [0,M1).

Inline propensity index I. Integer (3 X M1)
  I(:,i) gives the indices used by reaction #i. Where i is in [0,M1].

Inline propensity subdomain S. Integer sparse (<=Ncells X M1)
  S(:,i) lists all subdomains in which the corresponding inline propensity
  is off.

tspan. Double vector.
  Output times. tspan[0] is the start time and tspan[tlen-1] is the
  stop time.

vol. Double vector (length Ncells).
  vol[i] gives the volume of cell #i.

ldata. Double matrix (dsize X Ncells).
  Local data matrix, ldata(:,j) gives the local data vector for cell #j.

gdata. Double vector
  Global data matrix, passed to preopensity functions. 

sd. Integer vector (length Ncells).
  Subdomain number. sd[i] is the subdomain of cell #i. The vector sd
  can also be used to separate boundaries, line segments and points.

Format of sparse matrices:
  G, N and S are sparse matrices in compressed column format (CCS). D is sparse
  but in compressed row format (CRS), or equivalently, a transposed matrix in
  CCS format.
  jcD, irD, prD (double *)
  jcN, irN, prN (int *)
  jcG, irG (int *)

Propensities:
  a vector of function pointers (length Mreactions) is input by
  linking with the prototypes in propensities.h and function
  definitions in a user-specified .c-file. The type of this vector is
  PropensityFun which defines the input to a property function. See
  propensities.h for more details.

Ordering of the dofs:
  Dof #i is located in cell #(i/Mspecies), and the dofs located in
  cell #j is u0(:,j). Thus, u0 is understood as a matrix of size
  Mspecies X Ncells.

The output is a matrix U (Ndofs X length(tspan)).
  U(:,j) contains the state of the system at tspan(j).
*/
{
  double tt;
  double rdelta,rrdelta;
  double rand,cum,old;
  double *srrate,*rrate;
  double *sdrate,*Ddiag;
  double *rtimes;
  double old_rrate,old_drate;
  double totrate;

  int *node,*heap,*xx;
  int errcode = 0;
  long total_diffusion = 0,total_reactions = 0;
  int dof,col;
        
  int subvol,event,re,spec,to_spec;
  size_t i,j,it,k,timeix;
  size_t to_node,to_vol = 0;
  const size_t Ndofs = Ncells*Mspecies;

  ReportFun report = &URDMEreportFun;

  /* Main allocation. */

  /* current state vector */
  xx = (int *)MALLOC(Ndofs*sizeof(int));

  /* Reaction rate matrix (Mreactions X Ncells) and total rate
     vector. In rrate we store all propensities for chemical rections,
     and in srrate the sum of propensities in every subvolume. */
  rrate = (double *)MALLOC(Mreactions*Ncells*sizeof(double));
  srrate = (double *)MALLOC(Ncells*sizeof(double));

  /* Total diffusion rate vector (length Mcells). It will hold the
     total diffusion rates in each subvolume. */
  sdrate = (double *)MALLOC(Ncells*sizeof(double));

  /* The diagonal value of the D-matrix is used frequently. For
     efficiency, we store the negative of D's diagonal in Ddiag. */
  Ddiag = (double *)MALLOC(Ndofs*sizeof(double));
  for (i = 0; i < Ndofs; i++) {
    Ddiag[i] = 0.0;
    for (j = jcD[i]; j < jcD[i+1]; j++)
      if (irD[j] == i) {
	Ddiag[i] = -prD[j];
	break;
      }
  }

  /* Binary (min)heap. */
  rtimes = (double *)MALLOC(Ncells*sizeof(double));
  node = (int *)MALLOC(Ncells*sizeof(int));
  heap = (int *)MALLOC(Ncells*sizeof(int));

  /* Loop over Nreplicas cases. */
  for (k = 0; k < Nreplicas; k++) {
    it = 0;
    tt = tspan[0];
    old_rrate = 0.0;
    old_drate = 0.0;

    /* Set xx to the initial state. */
    memcpy(xx,&u0[k*Ndofs],Ndofs*sizeof(int));

    /* Pointer into data_time, ldata_time, gdata_time. */
    timeix = 0;

    /* set new master seed */
    srand48(seed_long[k]);

    /* Calculate the propensity for every reaction and every
       subvolume. Store the sum of the reaction intensities in each
       subvolume in srrate. */
    for (i = 0; i < Ncells; i++) {
      srrate[i] = 0.0;
      for (j = 0; j < M1; j++) {
	rrate[i*Mreactions+j] = 
	  inlineProp(&xx[i*Mspecies],&K[j*3],&I[j*3],&prS[jcS[j]],
		     jcS[j+1]-jcS[j],vol[i],sd[i]);
	srrate[i] += rrate[i*Mreactions+j];
      }
      for (; j < Mreactions; j++) {
	rrate[i*Mreactions+j] = 
	  (*rfun[j-M1])(&xx[i*Mspecies],tt,vol[i],&ldata[i*ldsize],gdata,
			&ldata_time[(timeix*Ncells+i)*ldtsize],
			&gdata_time[timeix*gdtsize],
			sd[i]);
	srrate[i] += rrate[i*Mreactions+j];
      }
    }

    /* Calculate the total diffusion rate for each subvolume. */
    for(i = 0; i < Ncells; i++) {
      sdrate[i] = 0.0;
      for(j = 0; j < Mspecies; j++)
	sdrate[i] += Ddiag[i*Mspecies+j]*xx[i*Mspecies+j];
    }

    /* Calculate times to next event (reaction or diffusion) in each
       subvolume and initialize heap. */
    for (i = 0; i < Ncells; i++) {
      rtimes[i] = -log(1.0-drand48())/(srrate[i]+sdrate[i])+tspan[0];
      heap[i] = node[i] = i;
    }
    initialize_heap(rtimes,node,heap,Ncells);

    /* Main loop. */
    for ( ; ; ) {
      /* Get the subvolume in which the next event occurred. This
	 subvolume is on top of the heap. */
      tt = rtimes[0];
      subvol = node[0];

      /* rather increase timeix first? */
      if (timeix+1 < dtlen && tt >= data_time[timeix+1]) {
	tt = data_time[++timeix];
 
	/* ldata_time and/or gdata_time was updated: recalculate
	   srrate and rrate in all voxels using the dependency
	   graph */
	for (subvol = 0; subvol < Ncells; subvol++) {
	  for (i = jcG[Mspecies+Mreactions], rdelta = 0.0; i < jcG[Mspecies+Mreactions+1]; i++) {
	    old = rrate[subvol*Mreactions+irG[i]];
	    j = irG[i];
	    if (j < M1)
	      rdelta += (rrate[subvol*Mreactions+j] = 
			 inlineProp(&xx[subvol*Mspecies],
				    &K[j*3],&I[j*3],&prS[jcS[j]],
				    jcS[j+1]-jcS[j],vol[subvol],sd[subvol]))-old;
	    else
	      rdelta += 
		(rrate[subvol*Mreactions+j] = 
		 (*rfun[j-M1])(&xx[subvol*Mspecies],tt,vol[subvol],
			       &ldata[subvol*ldsize],gdata,
			       &ldata_time[(timeix*Ncells+subvol)*ldtsize],
			       &gdata_time[timeix*gdtsize],
			       sd[subvol]))-old;
	  }
	  old = srrate[subvol]+sdrate[subvol];
	  srrate[subvol] += rdelta;
	  totrate = old+rdelta;
	  if (totrate > 0.0) {
	    if (!isinf(rtimes[heap[subvol]]))
	      rtimes[heap[subvol]] = old/totrate*(rtimes[heap[subvol]]-tt)+tt;
	    else
	      /* generate a new waiting time */
	      rtimes[heap[subvol]] = -log(1.0-drand48())/totrate+tt;
	  }
	  else
	    rtimes[heap[subvol]] = INFINITY;

	  update(heap[subvol],rtimes,node,heap,Ncells);
	}
	subvol = -1; /* signal nil event */
      }

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

      /* catch nil event: done reporting, so go for next event */
      if (subvol == -1) continue;

      /* First check if it is a reaction or a diffusion event. */
      totrate = srrate[subvol]+sdrate[subvol];
      rand = drand48();

      if (rand*totrate <= srrate[subvol]) {
	/* Reaction event. */
	event = -1;

	/* a) Determine the reaction re that did occur (direct SSA). */
	rand *= totrate;
	for (re = 0, cum = rrate[subvol*Mreactions]; 
	     re < Mreactions && rand > cum; 
	     re++, cum += rrate[subvol*Mreactions+re]);

	/* elaborate floating point fix: */
	if (re >= Mreactions) re = Mreactions-1;
	if (rrate[subvol*Mreactions+re] == 0.0) {
	  /* go backwards and try to find first nonzero reaction rate */
	  for ( ; re > 0 && rrate[subvol*Mreactions+re] == 0.0; re--);

	  /* No nonzero rate found, but a reaction was sampled. This can
	     happen due to floating point errors in the iterated
	     recalculated rates. */
	  if (rrate[subvol*Mreactions+re] == 0.0) {
	    /* nil event: zero out and move on */
	    srrate[subvol] = 0.0;
	    event = 0;
	    goto next_event;
	  }
	}

	/* b) Update the state of the subvolume subvol and sdrate[subvol]. */
	for (i = jcN[re]; i < jcN[re+1]; i++) {
	  xx[subvol*Mspecies+irN[i]] += prN[i];
	  if (xx[subvol*Mspecies+irN[i]] < 0) errcode = 1;
	  sdrate[subvol] += Ddiag[subvol*Mspecies+irN[i]]*prN[i];
	}

	/* c) Recalculate srrate[subvol] using dependency graph. */
	for (i = jcG[Mspecies+re], rdelta = 0.0; i < jcG[Mspecies+re+1]; i++) {
	  old = rrate[subvol*Mreactions+irG[i]];
	  j = irG[i];
	  if (j < M1)
	    rdelta += (rrate[subvol*Mreactions+j] = 
		       inlineProp(&xx[subvol*Mspecies],
				  &K[j*3],&I[j*3],&prS[jcS[j]],
				  jcS[j+1]-jcS[j],vol[subvol],sd[subvol]))-old;
	  else
	    rdelta += 
	      (rrate[subvol*Mreactions+j] = 
	       (*rfun[j-M1])(&xx[subvol*Mspecies],tt,vol[subvol],
			     &ldata[subvol*ldsize],gdata,
			     &ldata_time[(timeix*Ncells+subvol)*ldtsize],
			     &gdata_time[timeix*gdtsize],
			     sd[subvol]))-old;
	}
	srrate[subvol] += rdelta;

	total_reactions++; /* counter */
      }
      else {
	/* Diffusion event. */
	event = 1;

	/* a) Determine which species... */
	rand *= totrate;        
	rand -= srrate[subvol];
	for (spec = 0, dof = subvol*Mspecies, cum = Ddiag[dof]*xx[dof]; 
	     spec < Mspecies && rand > cum;
	     spec++, cum += Ddiag[dof+spec]*xx[dof+spec]);

	/* elaborate floating point fix: */
	if (spec >= Mspecies) spec = Mspecies-1;
	if (xx[dof+spec] == 0) {
	  /* go backwards and try to find first nonzero species */
	  for ( ; spec > 0 && xx[dof+spec] == 0; spec--);

	  /* No species to diffuse, but nonzero rate. This can happen
	     due to floating point errors in the iterated recalculated
	     diffusion rates. */
	  if (xx[dof+spec] == 0) {
	    /* nil event: zero out and move on */
	    sdrate[subvol] = 0.0;
	    event = 0;
	    goto next_event;
	  }
	}

	/* b) ...and then the direction of diffusion. */
	col = dof+spec;
	rand = drand48()*Ddiag[col];

	/* Search for diffusion direction. */
	for (i = jcD[col], cum = 0.0; i < jcD[col+1]; i++)
	  if (irD[i] != col && (cum += prD[i]) > rand) break;

	/* simple floating point fix: */
	if (i >= jcD[col+1]) i = jcD[col+1]-1;

	/* note: both pairs (subvol,to_vol) and (spec,to_spec) are
	   allowed to be non-equal */
	to_node = irD[i];
	to_vol = to_node/Mspecies;
	to_spec = to_node%Mspecies;
	/* hence this errorcode is no longer used: */
	/* if (subvol != to_vol && spec != to_spec) errcode = -1; */

	/* c) Execute the diffusion event (check for negative elements). */
	xx[subvol*Mspecies+spec]--;
	if (xx[subvol*Mspecies+spec] < 0) errcode = 2;
	xx[to_node]++;

	/* Save reaction and diffusion rates. */
	old_rrate = srrate[to_vol];
	old_drate = sdrate[to_vol];

	/* d.i) Recalculate the reaction rates in subvolume subvol
	   that are affected by species spec changing, using the
	   dependency graph G. */
	for (i = jcG[spec], rdelta = 0.0; i < jcG[spec+1]; i++) {
	  old = rrate[subvol*Mreactions+irG[i]];
	  j = irG[i];

	  if (j < M1)
	    rdelta += (rrate[subvol*Mreactions+j] = 
		       inlineProp(&xx[subvol*Mspecies],
				  &K[j*3],&I[j*3],
				  &prS[jcS[j]],jcS[j+1]-jcS[j],
				  vol[subvol],sd[subvol]))-old;
	  else
	    rdelta += 
	      (rrate[subvol*Mreactions+j] = 
	       (*rfun[j-M1])(&xx[subvol*Mspecies],tt,vol[subvol],
			     &ldata[subvol*ldsize],gdata,
			     &ldata_time[(timeix*Ncells+subvol)*ldtsize],
			     &gdata_time[timeix*gdtsize],
			     sd[subvol]))-old;
	}
	srrate[subvol] += rdelta;

	/* Adjust diffusion rate. */
	sdrate[subvol] -= Ddiag[subvol*Mspecies+spec];

	/* d.ii) Recalculate the reaction rates in subvolume to_vol
	   that are affected by species to_spec changing, using the
	   dependency graph G. */
	for (i = jcG[to_spec], rdelta = 0.0; i < jcG[to_spec+1]; i++) {
	  old = rrate[to_vol*Mreactions+irG[i]];
	  j = irG[i];

	  if (j < M1)
	    rdelta += (rrate[to_vol*Mreactions+j] = 
		       inlineProp(&xx[to_vol*Mspecies],
				  &K[j*3],&I[j*3],&prS[jcS[j]],
				  jcS[j+1]-jcS[j],
				  vol[to_vol],sd[to_vol]))-old;
	  else
	    rdelta += 
	      (rrate[to_vol*Mreactions+j] = 
	       (*rfun[j-M1])(&xx[to_vol*Mspecies],tt,vol[to_vol],
			     &ldata[to_vol*ldsize],gdata,
			     &ldata_time[(timeix*Ncells+to_vol)*ldtsize],
			     &gdata_time[timeix*gdtsize],
			     sd[to_vol]))-old;
	}
	srrate[to_vol] += rdelta;

	/* Adjust diffusion rate. */
	sdrate[to_vol] += Ddiag[to_node];

	total_diffusion++; /* counter */
      }

    next_event:
      /* Compute time to new event for this subvolume. */
      totrate = srrate[subvol]+sdrate[subvol];  
      if (totrate > 0.0)
	rtimes[0] = -log(1.0-drand48())/totrate+tt;
      else
	rtimes[0] = INFINITY;

      /* Update the heap. */
      update(0,rtimes,node,heap,Ncells);

      /* If it was a diffusion event (into another subvolume), also
	 update the other affected node. */
      if (event == 1 && subvol != to_vol) {
	totrate = srrate[to_vol]+sdrate[to_vol];      
	if (totrate > 0.0) {
	  if (!isinf(rtimes[heap[to_vol]]))
	    rtimes[heap[to_vol]] = 
	      (old_rrate+old_drate)/totrate*(rtimes[heap[to_vol]]-tt)+tt;
	  else
	    /* generate a new waiting time */
	    rtimes[heap[to_vol]] = -log(1.0-drand48())/totrate+tt;
	} 
	else
	  rtimes[heap[to_vol]] = INFINITY;

	update(heap[to_vol],rtimes,node,heap,Ncells);
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

  FREE(heap);
  FREE(node);
  FREE(rtimes);
  FREE(Ddiag);
  FREE(sdrate);
  FREE(srrate);
  FREE(rrate);
  FREE(xx);
}
/*----------------------------------------------------------------------*/
