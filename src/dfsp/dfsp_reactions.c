/* dfsp_reactions.c */

/* B. Drawert    2010-05-25. */
//#define _DEBUG

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "propensities.h"
#include "report.h"

int dfsp_reactions(int *xx, 
            const size_t *irN, const size_t *jcN, const int *prN,
            const size_t *irG, const size_t *jcG,
            PropensityFun *rfun,
            double tau_d,double tt,
            const double *vol, const double *data, const int *sd,
            const size_t Ncells,int Mspecies,int Mreactions,const size_t dsize)
/*

Input:

Ncells.
  Number of subvolumes.

Mspecies.
  Number of different species.

Mreactions.
  Total number of reactions.

dsize.
  Size of data vector sent to propensities.

Defines:
  Ndofs=Ncells*Mspecies.

Propensity prop. Vector of function pointers (length Mreactions).
  The type of this vector is PropensityFun which defines the input to a property
  function.
  typedef double (*PropensityFun)(const int *xx, double vol, const double
    *data, int sd);

Other inputs not covered above:
Ordering of the dofs: Dof #i is located in cell #(i/Mspecies),
and the dofs located in cell #j is u0(:,j). Thus, u0 is understood as a matrix
of size Mspecies X Ncells.

(cell -> vol) vol. Double vector (length Ncells).
  vol[i] gives the volume of cell #i.

(cell -> data) data. Double matrix (dsize X Ncells).
  Generalized data matrix, data(:,j) gives a data vector for cell #j.

(cell -> sd) sd. Integer vector (length Ncells).
  Subdomain number. sd[i] is the subdomain of cell #i. The vector sd can also
  be used to separate boundaries, line segments and points.

*/
{
    int i, j, k, re;

    double a0_rxn;
    double *a_rxn = (double *)malloc(Mreactions*sizeof(double));
    double tti,tau;
    double rand,cum,old,delta;

    int rxn_count=0;

    // Process each voxel
    for (i=0; i<Ncells; i++) {
        a0_rxn=0.0;
        for (j=0; j<Mreactions; j++) {
            a0_rxn += a_rxn[j] = (*rfun[j])(&xx[i*Mspecies],tt, vol[i], &data[i*dsize], sd[i]);
        
        }
        // Main Loop
        tti=0.0;
        while(tti<tau_d){  // has this voxel reached the end of the time step
            // find time to next reaction
            tau= -log( 1.0 - drand48() )/a0_rxn;
            if(tti+tau > tau_d) break;
            tti=tti+tau;
            // find which reaction fired
            cum=0.0;
            rand= drand48()*a0_rxn;
            re=-1;
            for (j=0; j<Mreactions; j++) {
                cum+=a_rxn[j];
                if(rand<=cum){ re=j; break; }
            }
            if(re==-1 || re>Mreactions){
                fprintf(stderr,"ERROR: could not find which reaction fired.\n");
                exit(0);
            }
            // process reaction
            for (k=jcN[re]; k<jcN[re+1]; k++) {
                xx[i*Mspecies+irN[k]] += prN[k];
                if (xx[i*Mspecies+irN[k]]<0){
                    fprintf(stderr,"ERROR: Negative species population, while processing reaction %i in voxel %i\n",re,i);
                    fprintf(stderr,"\trand=%e cum=%e a0_rxn=%e a_rxn[%i]=%e time=%e  tau=%e\n",rand,cum,a0_rxn,re,a_rxn[re],tti,tau);
                    for(j=0; j<Mreactions; j++) {
                        fprintf(stderr,"\t\ta_rxn[%i] = %e\n",j,a_rxn[j]);
                    }
                    exit(0);
                }
            }
            // update propensities using connectivity matrix
            for (k=jcG[Mspecies+re], delta=0.0; k<jcG[Mspecies+re+1]; k++) {
                j=irG[k];
                old = a_rxn[j];
                a_rxn[j]  =(*rfun[j])(&xx[i*Mspecies],tt+tti, vol[i], &data[i*dsize], sd[i]);
                delta = a_rxn[j]-old;
                a0_rxn += delta;
                if(a0_rxn < 0.0){
                    a0_rxn=0.0;
                }
            }
            // update rxn count
            rxn_count++;
        } //end Main Loop

    }
    free(a_rxn);
    return rxn_count;
    //--------------------------
	
}

