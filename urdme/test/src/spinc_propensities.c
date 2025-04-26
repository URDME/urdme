/* Propensity file for test "spin4". */

/* P. Bauer 2012-03-28 */
/* J. Cullhed 2008-08-04 */

#include <stdlib.h>
#include "propensities.h"

#define NR 8
#define N_TESTS 10

/* External variable used to select correct propensity function */
extern int testNumber;
        
/* Propensities. */
/*---------------------------------------------------------------------*/
double rFun1(const int *x, double t, double vol, const double *data, int sd)
/* This is the reaction X+Y+Z > xyz/vol^2 > @. * Not allowed in SSA! */
{
  return x[0]*x[1]*x[2]/(vol*vol);
}
/*---------------------------------------------------------------------*/
double rFun2(const int *x, double t, double vol, const double *data, int sd)
/* This is the reaction X+Y > xy/vol > @. */_
{
  return (double)x[0]*x[1]/vol;
}
/*---------------------------------------------------------------------*/
double rFun3(const int *x, double t, double vol, const double *data, int sd)
/* This is the reaction X+X > x(x-1) > @. */
{
  return (double)x[0]*(x[0]-1)/vol;
}
/*---------------------------------------------------------------------*/
double rFun4(const int *x, double t, double vol, const double *data, int sd)
/* This is the reaction Y+Y > y(y-1) > @. */
{
  return (double)x[1]*(x[1]-1)/vol;
}
/*---------------------------------------------------------------------*/
double rFun5(const int *x, double t, double vol, const double *data, int sd)
/* This is the reaction Z > z > @. */
{
  return (double)x[2];
}
/*---------------------------------------------------------------------*/
double rFun6(const int *x, double t, double vol, const double *data, int sd)
/* This is the reaction X+Y > xy > @. Only active in subdomain one. */
{
  if (sd==1) return (double)x[0]*x[1]/vol;
  return 0.0;
}
/*---------------------------------------------------------------------*/
double rFun7(const int *x, double t, double vol, const double *data, int sd)
/* This is the reaction Z > z > @. Only active in subdomain zero. */
{
  if (sd==0) return (double)x[2]/vol;
  return 0.0;
}
/*---------------------------------------------------------------------*/
double rFun8(const int *x, double t, double vol, const double *data, int sd)
/* Currently unused. */
{
  return 0.0;
}
/*---------------------------------------------------------------------*/

PropensityFun *ALLOC_propensities(void)
{
  PropensityFun *ptr = (PropensityFun *)malloc(sizeof(PropensityFun)*NR);
    
  switch(testNumber)
  {         
    case 1:
      ptr[0]=rFun2;
      break;
    case 2:
      ptr[0]=rFun3;
      ptr[1]=rFun4;
      break;
    case 3:
      ptr[0]=rFun3;
      ptr[1]=rFun4;
      ptr[2]=rFun5;
      break;
    case 4:
      ptr[0]=rFun6;
      ptr[1]=rFun7;
      break;
    case 5:
      ptr[0]=rFun1;
      break;
  }
  return ptr;
}


void FREE_propensities(PropensityFun* ptr)
{
  free(ptr);
}
