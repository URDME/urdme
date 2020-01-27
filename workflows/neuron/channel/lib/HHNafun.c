/* HHNafun.c - Auxiliary function definitions for H-H Na channels. */

/* S. Engblom 2019-12-03 */

#include <math.h>
#include "HHNafun.h"

double HHNa_alpha1(const double v)
{
  double t0 = 0.0;
  if (fabs(v-25) < 1e-2) {
    double p = v-25;
    
    t0 = v*5.0E-2-2.5E-1;
    p *= p;
    t0 += 8.333333333333333E-4*p;
    p *= p;
    t0 -= 1.388888888888889E-7*p;
  }
  else
    t0 = 0.1*(25-v)/(exp((25-v)/10)-1.0);

  return t0;
}

double HHNa_beta1(const double v)
{
  return 4.0*exp(-v/18);
}

double HHNa_alpha2(const double v)
{
  return 0.07*exp(-v/20);
}

double HHNa_beta2(const double v)
{
  return 1.0/(exp(3-0.1*v)+1.0);
}
