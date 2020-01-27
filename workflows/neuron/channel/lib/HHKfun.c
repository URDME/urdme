/* HHKfun.c - Auxiliary function definitions for H-H K channels. */

/* S. Engblom 2019-12-03 */

#include <math.h>
#include "HHKfun.h"

double HHK_alpha1(const double v)
{
  double t0 = 0.0;
  if (fabs(v-10) < 1e-2) {
    double p = v-10;
    t0 = 5.0E-2*(1+0.1*v);
    p *= p;
    t0 += 8.333333333333333E-5*p;
    p *= p;
    t0 -= 1.388888888888889E-8*p;
  }
  else
    t0 = 0.01*(10-v)/(exp(1-0.1*v)-1.0);

  return t0;
}

double HHK_beta1(const double v)
{
  return 0.125*exp(-v/80);
}
