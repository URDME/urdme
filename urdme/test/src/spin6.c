/* Propensity file for test "spin6". */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* S. Engblom 2017-02-18 (Revision) */
/*  P. Bauer 2012-03-28 */

#include "propensities.h"

#define NR 5

/* Reaction propensities. */
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  if (sd != 2)
    return x[3];
  else
    return 0.0;
}

double rFun2(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return x[0]*(x[0]-1)/vol;
}

double rFun3(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return x[1]*(x[1]-1)/vol;
}

double rFun4(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return x[2]*(x[2]-1)/vol;
}

double rFun5(const URDMEstate_t *x,double time,double vol,
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
  
  ptr[0] = rFun1;
  ptr[1] = rFun2;
  ptr[2] = rFun3;
  ptr[3] = rFun4;
  ptr[4] = rFun5;
    
  return ptr;
}

void FREE_propensities(PropensityFun* ptr)
{
  FREE(ptr);
}

void GET_jacobian(double *val,const URDMEstate_t *x,double time,double vol,
		  const double *ldata,const double *gdata,
		  const double *ldata_time,const double *gdata_time,int sd)
{
  val[0] = x[0]/vol+(x[0]-1.0)/vol;
  val[1] = x[1]*x[2]*1.0/(vol*vol);
  val[2] = x[1]/vol+(x[1]-1.0)/vol;
  val[3] = x[0]*x[2]*1.0/(vol*vol);
  val[4] = x[2]/vol+(x[2]-1.0)/vol;
  val[5] = x[0]*x[1]*1.0/(vol*vol);
  if (sd != 2)
    val[6] = 1.0;
  else
    val[6] = 0.0; 
}
