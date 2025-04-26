/* Propensity file for test "spin1". */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* S. Engblom 2017-02-16 */
/* P. Bauer 2012-03-26 */

#include "propensities.h"
#include "report.h"

/* number of reactions */
const int NR = 1;

/* forward declaration */
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd);

/* static propensity vector */
static PropensityFun ptr[] = {rFun1};

/*----------------------------------------------------------------------*/
/* Reaction propensities. */
double rFun1(const URDMEstate_t *x,double time,double vol,
             const double *ldata,const double *gdata,
             const double *ldata_time,const double *gdata_time,int sd)
{
  return x[0]*x[1]/vol;
}
/*----------------------------------------------------------------------*/
PropensityFun *ALLOC_propensities(size_t Mreactions)
{
  if (Mreactions > NR)
    PERROR("Wrong number of reactions.");
  return ptr;
}
/*----------------------------------------------------------------------*/
void FREE_propensities(PropensityFun* ptr)
{
  /* do nothing */
}
/*----------------------------------------------------------------------*/
