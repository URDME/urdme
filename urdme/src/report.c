/* Report function used with rdme_solve. */

/* A. Hellander 2008-12-15 */
/* B. Drawert 2010-12-12 */

// all files must compile without Matlab 
#include <stdlib.h>
#ifndef URDME_LIBMAT   
#include "mex.h"
#define printf mexPrintf
#define error_fn mexErrMsgTxt
#else //-------------------
#include <stdio.h>
#define error_fn(a) fprintf(stderr,a);exit(0);
#endif

void reportFun1(double time, const double t0, const double tend, long int total_diffusion,
		long int total_reaction, int errcode, int report_level)
/* Report function passed as argument to rdme_solve. */
{
  if (!errcode) {
    double prct=(time-t0)/(tend-t0)*100.0;
    
    switch (report_level) {
    case 1:
      printf("%i%% done.\n",(int)prct);
      break;
    case 2:
      printf("%i%% done.\n",(int)prct);
      printf("\t#Reaction  events = %li\n", total_reaction);
	  printf("\t#Diffusion events = %li\n", total_diffusion);
	  break;
    default:
      /* May be expanded in the future... */
      break;
    }
  } else {
    switch (errcode) {
    case 1:
      error_fn("Negative state detected.");
      break;
    default:
      error_fn("Unknown error code.");
      break;
    }
  }
}
