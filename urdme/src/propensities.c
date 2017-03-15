/* propensities.c - Empty declaration of URDME propensity functions. */

/* S. Engblom 2017-02-17 (Major revision, URDME 1.3, Comsol 5) */

#include "mex.h"
#include "propensities.h"

/* This file should be compiled and linked in case no propensity file
   is used. */

/* placeholders */
PropensityFun *ALLOC_propensities(size_t Mreactions)
{
  if (Mreactions > 0)
    mexErrMsgTxt("Wrong number of reactions.");
  return NULL;
}
void FREE_propensities(PropensityFun* ptr) { /* do nothing */ }
