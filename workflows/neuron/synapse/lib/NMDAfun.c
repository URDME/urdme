/* NMDAfun.c - Auxiliary function definitions for NMDA synaptic propensities. */

/* S. Engblom 2019-11-29 */

#include <math.h>
#include "NMDAfun.h"

double rmb(const double v)
{
  return exp((v-40)*valence*memb_fraction/25);
}

double rmu(const double v)
{
  return exp(-(v-40)*valence*(1-memb_fraction)/25);
} 
