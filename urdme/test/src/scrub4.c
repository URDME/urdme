/* Propensity file for test "scrub4". */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* S. Engblom 2017-02-22 */
/* P. Bauer 2012-03-28 */

#include "propensities.h"

/* number of reactions */
const int NR = 3;

/* forward declaration */
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd);
double rFun2(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd);
double rFun3(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd);

/* static propensity vector */
static PropensityFun ptr[] = {rFun1,rFun2,rFun3};

/*----------------------------------------------------------------------*/
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  if (sd == 3)
    return 1000.0*vol;
  else
    return 0.0;
}
/*----------------------------------------------------------------------*/
double rFun2(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  if (sd == 4)
    return (double)x[0];
  else
    return 0.0;
}
/*----------------------------------------------------------------------*/
double rFun3(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  if (sd == 2)
    return (double)x[1];
  else
    return 0.0;
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
