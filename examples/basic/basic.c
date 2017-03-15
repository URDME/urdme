/* Propensity file for basic URDME example.
 *
 * Reaction X+Y <--> Z in a 3D sphere.
 *
 * S. Engblom 2017-02-19 (Revision)
 * S. Engblom 2013-03-27 (Minor revision)
 * P. Bauer 2013-03-20
 */

#include "propensities.h"

/* number of reactions */
const int NR = 2;

/* forward declarations */
double rFun1(const int *x, double t, double vol,
	     const double *ldata, const double *gdata,int sd);
double rFun2(const int *x, double t, double vol,
	     const double *ldata, const double *gdata,int sd);

/* static propensity vector */
static PropensityFun ptr[] = {rFun1,rFun2};

/*----------------------------------------------------------------------*/
double rFun1(const int *x, double t, double vol,
	     const double *ldata, const double *gdata, int sd)
/* X+Y --> Z */
{
  return x[0]*x[1]/vol; 
}
/*----------------------------------------------------------------------*/
double rFun2(const int *x, double t, double vol,
	     const double *ldata, const double *gdata, int sd)
/* Z --> X+Y */
{
  return x[2];
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
