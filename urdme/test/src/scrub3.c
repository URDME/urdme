/* Propensity file for test "scrub3". */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* S. Engblom 2017-02-22 */
/* P. Bauer 2012-03-28 */

#include "propensities.h"

/* number of reactions */
const int NR = 4;

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
double rFun4(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd);

/* static propensity vector */
static PropensityFun ptr[] = {rFun1,rFun2,rFun3,rFun4};

/*----------------------------------------------------------------------*/
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  if (sd!=1)
    return 10000*vol;
  return 0.0;
}
/*----------------------------------------------------------------------*/
double rFun2(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  if (sd!=1)
    return 10000*vol;
  return 0.0;
}
/*----------------------------------------------------------------------*/
double rFun3(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return 0.01*x[0]*x[1]/vol;
}
/*----------------------------------------------------------------------*/
double rFun4(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return 0.01*x[2]*x[4]/vol;
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
