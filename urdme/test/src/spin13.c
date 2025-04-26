/* Propensity file for test "spin11". */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* P. Bauer 2013-04-10 */

#include "propensities.h"

#define NR 4

/* Reaction propensities. */
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return 2*x[0];
}

double rFun2(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return 5*x[1];
}

double rFun3(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return 8*x[2];
}

double rFun4(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return 17*x[3];
}

PropensityFun *ALLOC_propensities(size_t Mreactions)
{
  if (Mreactions > NR)
    mexErrMsgTxt("Wrong number of reactions.");
  PropensityFun *ptr = MALLOC(sizeof(PropensityFun)*NR);
  
  ptr[0] = rFun1;
  ptr[1] = rFun2;
  ptr[2] = rFun3;
  ptr[3] = rFun4;
    
  return ptr;
}

void FREE_propensities(PropensityFun* ptr)
{
  FREE(ptr);
}
