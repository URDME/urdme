/* Propensity file for test "scrub1". */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* S. Engblom 2017-02-22 */
/* P. Bauer 2012-03-28 */

#include "propensities.h"

/* number of reactions */
const int NR = 1;

/* forward declaration */
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd);

/* static propensity vector */
static PropensityFun ptr[] = {rFun1};

/*----------------------------------------------------------------------*/
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return 1.9293e-14*x[0]*x[1]/vol;
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
