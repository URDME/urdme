/* report.c - URDME report function. */

/* S. Engblom 2017-02-17 (Major revision, URDME 1.3, Comsol 5) */
/* S. Engblom 2014-06-10 (Revision) */
/* B. Drawert 2010-12-12 */
/* A. Hellander 2008-12-15 */

#include "mex.h"
#include "report.h"

/*----------------------------------------------------------------------*/
void URDMEreportFun(double time,double t0,double tend,
		    long total_diffusion,long total_reaction,
		    int errcode,int report_level)
{
  static int ncalls;

  if (!errcode) {
    double prct = (time-t0)/(tend-t0)*100.0;
    
    switch (report_level) {
    case 1:
      PRINTF("%i%% ",(int)prct);
      if (++ncalls == 20) {
	PRINTF("\n");
	ncalls = 0;
      }
      break;
    case 2:
      PRINTF("%i%% done.\n",(int)prct);
      PRINTF("\t#Reaction  events = %li\n",total_reaction);
      PRINTF("\t#Diffusion events = %li\n",total_diffusion);
      break;
    case 3:
      PRINTF("%i%% done.\n",(int)prct);
      PRINTF("\t#Reaction  events = %li\n",total_reaction);
      PRINTF("\t#Diffusion events = %li\n",total_diffusion);
      /* callback allows for interactive exit: */
      PRINTF("Type DBCONT (or RETURN) to continue, DBQUIT to quit.\n");
      mexCallMATLAB(0,NULL,0,NULL,"keyboard");
      break;
    default:
      /* May be expanded in the future... */
      break;
    }
  }
  else {
    switch (errcode) {
    case 1:
      PERROR("Negative state detected (reaction).\n");
      break;
    case 2:
      PERROR("Negative state detected (diffusion).\n");
      break;
    case -1:
      PERROR("A transport event can not both change voxel and species "
	     "at the same time.\n");
      break;
    default:
      PRINTF("Unknown error code = %d.\n",errcode);
      PERROR("Bailing out.\n");
      break;
    }
  }
}
/*----------------------------------------------------------------------*/
