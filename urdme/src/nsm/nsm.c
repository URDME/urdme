/* nsm.c - URDME NSM solver. */

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
	 const double *tspan,const size_t tlen,
	 int *U,
	 const double *vol,const double *ldata,const double *gdata,
	 const int *sd,
	 const size_t Ncells,
	 const size_t Mspecies,const size_t Mreactions,
	 const size_t dsize,int report_level,
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

report_level
  The desired degree of feedback during simulations. 0, 1, and 2 are
  currently supported options.

Diffusion matrix D. Double sparse (Ndofs X Ndofs).
  Macroscopic diffusion matrix. D(i,j) is the diffusion rate from dof #j to 
  dof #i. This matrix uses the CSR-format and not CSC because fast access to
  rows is needed.

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
  double tt = tspan[0];
  double rdelta,rrdelta;
  double rand,cum,old;
  double *srrate,*rrate;
  double *sdrate,*Ddiag;
  double *rtimes;
  double old_rrate = 0.0,old_drate = 0.0;
  double totrate;

  int *node,*heap,*xx;
  long int total_reactions = 0;
  long int total_diffusion = 0;
  int dof,col;
        
  int subvol,event,re,spec,errcode = 0;
  size_t i,j,it = 0;
  size_t to_node,to_vol = 0;
  const size_t Ndofs = Ncells*Mspecies;

  ReportFun report = &URDMEreportFun;

  /* Set xx to the initial state. */
  xx = (int *)MALLOC(Ndofs*sizeof(int));
  memcpy(xx,u0,Ndofs*sizeof(int));

  /* Create reaction rate matrix (Mreactions X Ncells) and total rate
     vector. In rrate we store all propensities for chemical rections,
     and in srrate the sum of propensities in every subvolume. */
  rrate = (double *)MALLOC(Mreactions*Ncells*sizeof(double));
  srrate = (double *)MALLOC(Ncells*sizeof(double));

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
	(*rfun[j-M1])(&xx[i*Mspecies],tt,vol[i],&ldata[i*dsize],gdata,sd[i]);
      srrate[i] += rrate[i*Mreactions+j];
    }
  }

  /* Total diffusion rate vector (length Mcells). It will hold the
     total diffusion rates in each subvolume. */
  sdrate = (double *)MALLOC(Ncells*sizeof(double));

  /* The diagonal value of the D-matrix is used frequently. For
     efficiency, we store the negative of D's diagonal in Ddiag. */
  Ddiag = (double *)MALLOC(Ndofs*sizeof(double));
  for (i = 0; i < Ndofs; i++) {
    Ddiag[i] = 0.0;
    for (j = jcD[i]; j < jcD[i+1]; j++)
      if (irD[j] == i) Ddiag[i] = -prD[j];
  }

  /* Calculate the total diffusion rate for each subvolume. */
  for(i = 0; i < Ncells; i++) {
    sdrate[i] = 0.0;
    for(j = 0; j < Mspecies; j++)
      sdrate[i] += Ddiag[i*Mspecies+j]*xx[i*Mspecies+j];
  }

  /* Create binary (min)heap. */
  rtimes = (double *)MALLOC(Ncells*sizeof(double));
  node = (int *)MALLOC(Ncells*sizeof(int));
  heap = (int *)MALLOC(Ncells*sizeof(int));

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

    /* Store solution if the global time counter tt has passed the
       next time in tspan. */
    if (tt >= tspan[it] || isinf(tt)) {
      for (; it < tlen && (tt >= tspan[it] || isinf(tt)); it++) {
	if (report_level)
	  report(tspan[it],tspan[0],tspan[tlen-1],
		 total_diffusion,total_reactions,0,report_level);
	memcpy(&U[Ndofs*it],xx,Ndofs*sizeof(int));
      }

      /* If the simulation has reached the final time, exit. */     
      if (it >= tlen) break;
    }

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
			   &ldata[subvol*dsize],gdata,sd[subvol]))-old;
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

      /* b) and then the direction of diffusion. */
      col = dof+spec;
      rand = drand48()*Ddiag[col];

      /* Search for diffusion direction. */
      for (i = jcD[col], cum = 0.0; i < jcD[col+1]; i++)
        if (irD[i] != col && (cum += prD[i]) > rand) break;

      /* simple floating point fix: */
      if (i >= jcD[col+1]) i = jcD[col+1]-1;

      /* note: subvol and to_vol are allowed to be equal */
      to_node = irD[i];
      to_vol = to_node/Mspecies;

      /* c) Execute the diffusion event (check for negative elements). */
      xx[subvol*Mspecies+spec]--;
      if (xx[subvol*Mspecies+spec] < 0) errcode = 2;
      xx[to_node]++;

      if (subvol != to_vol) {
	/* Save reaction and diffusion rates. */
	old_rrate = srrate[to_vol];
	old_drate = sdrate[to_vol];

	/* d) Recalculate the reaction rates using dependency graph G. */
	for (i = jcG[spec], rdelta = 0.0, rrdelta = 0.0; i < jcG[spec+1]; i++) {
	  old = rrate[subvol*Mreactions+irG[i]];
	  j = irG[i];

	  if (j < M1) {
	    rdelta += (rrate[subvol*Mreactions+j] = 
		       inlineProp(&xx[subvol*Mspecies],
				  &K[j*3],&I[j*3],
				  &prS[jcS[j]],jcS[j+1]-jcS[j],
				  vol[subvol],sd[subvol]))-old;
	    old = rrate[to_vol*Mreactions+j];
	    rrdelta += (rrate[to_vol*Mreactions+j] = 
			inlineProp(&xx[to_vol*Mspecies],
				   &K[j*3],&I[j*3],&prS[jcS[j]],
				   jcS[j+1]-jcS[j],
				   vol[to_vol],sd[to_vol]))-old;
	  }
	  else{
	    rdelta += 
	      (rrate[subvol*Mreactions+j] = 
	       (*rfun[j-M1])(&xx[subvol*Mspecies],tt,vol[subvol],
			     &ldata[subvol*dsize],gdata,sd[subvol]))-old;
	    old = rrate[to_vol*Mreactions+j];
	    rrdelta += 
	      (rrate[to_vol*Mreactions+j] = 
	       (*rfun[j-M1])(&xx[to_vol*Mspecies],tt,vol[to_vol],
			     &ldata[to_vol*dsize],gdata,sd[to_vol]))-old;
	  }
	}

	srrate[subvol] += rdelta;
	srrate[to_vol] += rrdelta;

	/* Adjust diffusion rates. */
	sdrate[subvol] -= Ddiag[subvol*Mspecies+spec];
	sdrate[to_vol] += Ddiag[to_node];
      }
      else {
	/* d) Recalculate the reaction rates using dependency graph G. */
	for (i = jcG[spec], rdelta = 0.0, rrdelta = 0.0; i < jcG[spec+1]; i++) {
	  old = rrate[subvol*Mreactions+irG[i]];
	  j = irG[i];

	  if (j < M1)
	    rrdelta += (rrate[subvol*Mreactions+j] = 
			inlineProp(&xx[subvol*Mspecies],
				   &K[j*3],&I[j*3],&prS[jcS[j]],
				   jcS[j+1]-jcS[j],
				   vol[subvol],sd[subvol]))-old;
	  else
	    rdelta += 
	      (rrate[subvol*Mreactions+j] = 
	       (*rfun[j-M1])(&xx[subvol*Mspecies],tt,vol[subvol],
			     &ldata[subvol*dsize],gdata,sd[subvol]))-old;
	}
	/*** formally a bug: it is allowed to diffuse from one
	     subvolume as "species 1" and enter another subvolume as
	     "species 2"; hence in the case subvol != to_vol, one
	     should rather loop over the union of G(:,spec) and
	     G(:,to_spec), with to_spec = to_node%Mspecies, at least
	     whenever to_spec != spec */

	srrate[subvol] += rdelta;

	/* Adjust diffusion rates. */
	sdrate[subvol] -= Ddiag[subvol*Mspecies+spec];
	sdrate[to_vol] += Ddiag[to_node];
      }

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
        
    /* Check for error codes. */
    if (errcode) {
      /* Report the error that occurred and exit. */
      memcpy(&U[Ndofs*it],xx,Ndofs*sizeof(int));
      report(tt,tspan[0],tspan[tlen-1],
	     total_diffusion,total_reactions,errcode,
	     report_level);
      break;
    }
  }

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
