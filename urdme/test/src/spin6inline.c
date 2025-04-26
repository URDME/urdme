/* Propensity file for test "spin6". */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* S. Engblom 2017-02-18 (Revision) */
/* P. Bauer 2012-03-28 */

#include "propensities.h"

#define NR 1

/* Reaction propensities. Only the cubic reaction remains here since
   the rest is made inline. */
double rFun0(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return x[0]*x[1]*x[2]/(vol*vol);
}

PropensityFun *ALLOC_propensities(size_t Mreactions)
{
  if (Mreactions > NR)
    mexErrMsgTxt("Wrong number of reactions.");
  PropensityFun *ptr = (PropensityFun *)MALLOC(sizeof(PropensityFun)*NR);
  
  ptr[0] = rFun0;
    
  return ptr;
}

void FREE_propensities(PropensityFun* ptr)
{
  FREE(ptr);
}

void GET_jacobian(double *val,
                  const URDMEstate_t *xstate,double time,double vol,
                  const double *ldata,const double *gdata,
                  const double *ldata_time,const double *gdata_time,int sd)
{
  val[0] = xstate[1]*xstate[2]*1.0/(vol*vol);
  val[1] = xstate[0]*xstate[2]*1.0/(vol*vol);
  val[2] = xstate[0]*xstate[1]*1.0/(vol*vol);
}
