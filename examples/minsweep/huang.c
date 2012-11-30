/* huang.c

   Model file for the MinD/MinE-sweep example. 
   The model is taken from Huang et al. (2003).

   S. Engblom 2011-06-01 (Revision)
   A. Hellander 2010-02-19. 
*/
#include <stdlib.h>
#include "propensities.h"

// ordering of species
enum {MinD_c_atp,MinD_m,MinD_e,MinDE,MinD_c_adp};

// domain indicator values
enum {CYTOSOL = 1,MEMBRANE};

// number of reactions
const int NR = 6;

// rate constants
const double NA  = 6.0221415e23;
const double sigma_d  = 2.5e-8;
const double sigma_dD = 0.0016e-18;
const double sigma_e  = 0.093e-18;
const double sigma_de = 0.7;
const double sigma_dt = 1.0;

// declaration of exported propensities
double rFun1(const int *,double,double,const double *,int);
double rFun2(const int *,double,double,const double *,int);
double rFun3(const int *,double,double,const double *,int);
double rFun4(const int *,double,double,const double *,int);
double rFun5(const int *,double,double,const double *,int);
double rFun6(const int *,double,double,const double *,int);

/*----------------------------------------------------------------------*/
PropensityFun *ALLOC_propensities(void)
/* Exported allocation. */
{
  PropensityFun *ptr = (PropensityFun *)malloc(sizeof(PropensityFun)*NR);

  ptr[0] = rFun1;
  ptr[1] = rFun2;
  ptr[2] = rFun3;
  ptr[3] = rFun4;
  ptr[4] = rFun5;
  ptr[5] = rFun6;

  return ptr;
}
/*----------------------------------------------------------------------*/
void FREE_propensities(PropensityFun* ptr)
/* Exported deallocation. */
{
  free(ptr);
}
/*----------------------------------------------------------------------*/

// now follows the definition of the reaction propensities

/*----------------------------------------------------------------------*/
double rFun1(const int *x,double t,double vol,const double *data,int sd)
/* MinD_c_atp -> MinD_m */
{
  if (sd == MEMBRANE)
    return sigma_d*x[MinD_c_atp]/data[0];
  return 0.0;
}
/*----------------------------------------------------------------------*/
double rFun2(const int *x,double t,double vol,const double *data,int sd)
/* MinD_c_atp + MinD_m -> 2MinD_m (cooperative binding) */
{
  if (sd == MEMBRANE)
    return sigma_dD*x[MinD_c_atp]*x[MinD_m]/vol;
  return 0.0;
}
/*----------------------------------------------------------------------*/
double rFun3(const int *x,double t,double vol,const double *data,int sd)
/* MinD_m + MinD_e -> MinDE */
{
  if (sd == MEMBRANE)
    return sigma_e*x[MinD_m]*x[MinD_e]/vol;
  return 0.0;
}
/*----------------------------------------------------------------------*/
double rFun4(const int *x,double t,double vol,const double *data,int sd)
/* MinDE -> MinD_c_adp + MinD_e */
{
  if (sd == MEMBRANE)
    return sigma_de*x[MinDE];
  return 0.0;
}
/*----------------------------------------------------------------------*/
double rFun5(const int *x,double t,double vol,const double *data,int sd)
/* MinD_c_adp -> MinD_c_atp */
{
  return sigma_dt*x[MinD_c_adp];
}
/*----------------------------------------------------------------------*/
double rFun6(const int *x,double t,double vol,const double *data,int sd)
/* MinDE + MinD_c_atp -> MinD_m + MinDE */
{
  if (sd == MEMBRANE)
    return sigma_dD*x[MinDE]*x[MinD_c_atp]/vol;
  return 0.0;
}
/*----------------------------------------------------------------------*/
