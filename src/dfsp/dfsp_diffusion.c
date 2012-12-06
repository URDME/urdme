
/* A. Hellander 2010-04-02.
   B. Drawert   2010-05-25 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "propensities.h"
#include "dfsp.h"

/* One timestep of dfsp diffusion, moving every molecule individually. */
int dfsp_diffusion(int *xx, const size_t *irD, const size_t *jcD, const double *prD,
                   const size_t Ncells,const int Mspecies)
{
	
	int i,j,dof,to_dof,nx,Ndofs; 
	double cumsum,rand;
	
	Ndofs=Ncells*Mspecies;
	
	/* Need copy of state vector. */
	int *xtemp;
	xtemp = (int *)malloc(Ndofs*sizeof(int));
	memcpy(xtemp,xx,Ndofs*sizeof(int));
	
	/* For every dof... */
	for (dof=0;dof<Ndofs;dof++)
	{
		/* Move one molecule at a time... */
		nx = xx[dof];
		
		for (j=0;j<nx;j++){
			/* throw a random number */		
			rand=drand48();	
			/* Sample diffusion event */
			for(i=jcD[dof],cumsum=0.0;i<jcD[dof+1];i++)
				if ((cumsum+=prD[i])>rand){
					break;
				}
			/* This shouldn't happen if columns sum exactly to 1. However, roundoff 
			   may cause that condition to be violated. */ 
			if (i>=jcD[dof+1]) i--;
			
			to_dof = irD[i];
			
			/* Update state */
			xtemp[dof]--;
			xtemp[to_dof]++;
		}
	}
		
	memcpy(xx,xtemp,Ndofs*sizeof(int));
	free(xtemp);
    return 1;
		
}
