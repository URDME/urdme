/* Propensity file for test "spin15". */
/* Linear birth-death model where parameters are passed in gdata. */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* S. Engblom 2017-02-25 */

#include "propensities.h"

/* number of reactions */
const int NR = 2;

/* forward declaration */
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
	     const double *ldata_time,const double *gdata_time,int sd);
double rFun2(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd);

/* static propensity vector */
static PropensityFun ptr[] = {rFun1,rFun2};

/*----------------------------------------------------------------------*/
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return gdata[0]*vol;
}
/*----------------------------------------------------------------------*/
double rFun2(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return gdata[1]*x[0];
}
/*----------------------------------------------------------------------*/
PropensityFun *ALLOC_propensities(size_t Mreactions)
{
  if (Mreactions > NR)
    mexErrMsgTxt("Wrong number of reactions.");

  return ptr;
}
/*----------------------------------------------------------------------*/
void FREE_propensities(PropensityFun* ptr)
{
  /* do nothing */
}
/*----------------------------------------------------------------------*/
