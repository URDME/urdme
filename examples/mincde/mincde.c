/* Model file for the mincde example. 
*
*  A. Hellander 2010-02-19. 
*
*  Model from Fange and Elf, 2006. 
* 
*/

#include <stdlib.h>
#include "propensities.h"

#define MinD_c_atp 0
#define MinD_m	   1
#define MinD_e	   2
#define MinDE	   3
#define MinD_c_adp 4

#define CYTOSOL  1
#define MEMBRANE 2

#define NR		 5

/* Rate constants. */

const double NA   = 6.022e23;
const double kd   = 1.25e-8;
const double kdd  = 9.0e6; 
const double kde  = 5.56e7;
const double ke   = 0.7;
const double k_adp= 1.0;

 
/* Reaction propensities. */
double rFun1(const int *x, double t, double vol, const double *data, int sd)
/* MinD_c_atp -> MinD_m */
{
	if (sd == MEMBRANE) 
		return kd*x[MinD_c_atp]/data[0]; 
	return 0.0;
}

double rFun2(const int *x, double t, double vol, const double *data, int sd)
/* MinD_c_atp + MinD_m -> 2MinD_m */
{
	if (sd == MEMBRANE) 
		return kdd*x[MinD_c_atp]*x[MinD_m]/(1000.0*NA*vol);
	return 0.0;
}

double rFun3(const int *x, double t, double vol, const double *data, int sd)
/* MinD_m + MinD_e -> MinDE */
{
   if (sd == MEMBRANE) 
	   return kde*x[MinD_m]*x[MinD_e]/(1000.0*NA*vol);
   return 0.0;
}

double rFun4(const int *x, double t, double vol, const double *data, int sd)
/* MinDE -> MinD_e + MinD_c_adp*/
{
  if (sd == MEMBRANE) 
	  return ke*x[MinDE];
  return 0.0;
}

double rFun5(const int *x, double t, double vol, const double *data, int sd)
/* MinD_c_adp -> MinD_c_atp */
{
   return k_adp*x[MinD_c_adp];
}


PropensityFun *ALLOC_propensities(void)
{
  PropensityFun *ptr = malloc(sizeof(PropensityFun)*NR);
  
  ptr[0]=rFun1;
  ptr[1]=rFun2;
  ptr[2]=rFun3;
  ptr[3]=rFun4;
  ptr[4]=rFun5;
    
  return ptr;
}

void FREE_propensities(PropensityFun* ptr)
{
  free(ptr);
}

