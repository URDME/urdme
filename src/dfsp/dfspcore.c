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
    xx = (int *)malloc(Ndofs*sizeof(int));
    memcpy(xx,u0,Ndofs*sizeof(int));
	
    /* Number of non-zeros entries in diffusion matrix. */
    nnz = jcD[Ndofs];
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
                report(tspan[it], tspan[0], tspan[tlen-1], total_diffusion, total_reactions, 0,report_level);
              }
              memcpy(&U[Ndofs*it], &xx[0], Ndofs*sizeof(int));
          }

          /* If the simulation has reached the final time, exit. */
          if (it>=tlen) break;

        }
		
        /* First order splitting... */ 
 
        /* do dfsp diffusion (one full timestep),  */
        total_diffusion += dfsp_diffusion(xx,irD,jcD,prD,Ncells,Mspecies);

        /* and then reactions (one full timestep) starting with xx. */
	    total_reactions += dfsp_reactions(xx,
            (const size_t *)irN, (const size_t *)jcN, prN,
            (const size_t *)irG, (const size_t *)jcG,
            rfun, tau_d, tt, vol, data, sd, Ncells, Mspecies, Mreactions, dsize);
        tt+=tau_d;
		
    }
	
    FREE_propensities(rfun);

    free(xx);
 
}

