/* Propensity file for test "spin8". */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* P. Bauer 2012-03-28 */

#include "propensities.h"

#define NR 3

/* Reaction propensities. */
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return x[0]*x[1]/vol;
}

double rFun2(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return 100*x[1]*x[2]/vol;
}

double rFun3(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return x[0]*x[2]/vol;
}

PropensityFun *ALLOC_propensities(size_t Mreactions)
{
  if (Mreactions > NR)
    mexErrMsgTxt("Wrong number of reactions.");
  PropensityFun *ptr = (PropensityFun *)MALLOC(sizeof(PropensityFun)*NR);
  
  ptr[0] = rFun1;
  ptr[1] = rFun2;
  ptr[2] = rFun3;
    
  return ptr;
}

void FREE_propensities(PropensityFun* ptr)
{
  FREE(ptr);
}
