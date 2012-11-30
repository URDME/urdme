/* dfsp.c */
/* A. Hellander  2010-04-05. */
/* B. Drawert    2010-05-25. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "propensities.h"
#include "report.h"
#include "dfsp.h"
#include "mex.h"

#if !defined(MALLOC) || !defined(FREE)
#error "Must define MALLOC and FREE."
#endif



void *MALLOC(size_t);
void FREE(void *);

void dfsp_core(const size_t *irD, const size_t *jcD, const double *prD,
            const int *u0, const size_t *irN, const size_t *jcN, const int *prN,
            const size_t *irG, const size_t *jcG, const double *tspan, const size_t tlen, int *U,
            const double *vol, const double *data, const int *sd, const size_t Ncells,
            const size_t Mspecies, const size_t Mreactions, const size_t dsize,
            int report_level, const double tau_d)
{


	int i,nnz;
    long int total_reactions=0;
    long int total_diffusion=0;

	double tt,tend;
	tt =   tspan[0];
	tend = tspan[tlen-1];

	int Ndofs = Ncells*Mspecies;

	PropensityFun *rfun;
	rfun =  ALLOC_propensities();

    int *xx;
    xx = MALLOC(Ndofs*sizeof(int));
    memcpy(xx,u0,Ndofs*sizeof(int));
	
    /* Number of non-zeros entries in diffusion matrix. */
    nnz = jcD[Ndofs];

    /* Make a copy of diffusion matrix with order optimized columns (for speed).
       The entries in the column is sorted in decending order, so that the most
       likely event if placed first in the array that is searched for the 
       event when sampling in the diffusion step. */
	/*double *s;
	size_t *ir,*jc;
	s  = MALLOC(nnz*sizeof(double));
	ir = MALLOC(nnz*sizeof(size_t));
	jc = MALLOC((Ndofs+1)*sizeof(size_t));
	memcpy(jc,jcD,(Ndofs+1)*sizeof(size_t));
	memcpy(ir,irD,nnz*sizeof(size_t));
	memcpy(s,prD,nnz*sizeof(double));

	int start=0,stop=0;
	for (i=0;i<Ndofs;i++){
		start = jcD[i];
		stop  = jcD[i+1];		
		insertion_sort(&s[start],&ir[start],stop-start);
	}*/

	int it=0;

    ReportFun report;
    if (report_level)
       report=&reportFun1;
    else
        report=NULL;

    for (;;)
    {
        /* Store solution if the global time counter tt has passed the next time is tspan. */
        if (tt>=tspan[it]||isinf(tt)) {
          for (; it<tlen && (tt>=tspan[it]||isinf(tt)); it++){
              if (report){
                //fprintf(stderr,"Calling report()\n");
                report(tspan[it], tspan[0], tspan[tlen-1], total_diffusion, total_reactions, 0,report_level);
              }
              //fprintf(stderr,"Saving state vector\n");
              memcpy(&U[Ndofs*it], &xx[0], Ndofs*sizeof(int));
          }

          /* If the simulation has reached the final time, exit. */
          if (it>=tlen) break;

        }
		
        /* First order splitting... */ 
 
        /* do dfsp diffusion (one full timestep),  */
        //fprintf(stderr,"Diffusion Step\n");
        total_diffusion += dfsp_diffusion(xx,irD,jcD,prD,Ncells,Mspecies);

        /* and then reactions (one full timestep) starting with xx. */
        //fprintf(stderr,"Reaction Step\n");
	    total_reactions += dfsp_reactions(xx, 
            (const size_t *)irN, (const size_t *)jcN, prN,
            (const size_t *)irG, (const size_t *)jcG,
            rfun, tau_d, tt, vol, data, sd, Ncells, Mspecies, Mreactions, dsize);	

        tt+=tau_d;
        //fprintf(stderr,"One full step completed, t=%e\n",tt);
		
    }
	
    FREE_propensities(rfun);

    //FREE(ir);
    //FREE(jc);
    //FREE(s);
    FREE(xx);
 
}

