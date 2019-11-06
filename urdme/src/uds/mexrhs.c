/* mexrhs.c - Mex-interface for use with the UDS solver in URDME. */

/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* S. Engblom 2017-02-24 */

#include <stdlib.h>

#include "mex.h"
#include "matrix.h"

#include "inline.h"
#include "propensities.h"

/*----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /* check syntax */
  if (nrhs != 10) mexErrMsgTxt("Wrong number of arguments.");

  /* load arguments */
  const double tt = mxGetScalar(prhs[0]);
  const mxArray *mxU0 = prhs[1];
  const size_t Mreactions = (size_t)mxGetScalar(prhs[2]);
  const mxArray *mxVol = prhs[3];
  const mxArray *mxLData = prhs[4];
  const mxArray *mxGData = prhs[5];
  const mxArray *mxSd = prhs[6];
  const mxArray *mxK = prhs[7];
  const mxArray *mxI = prhs[8];
  const mxArray *mxS = prhs[9];

  /* get problem dimensions */
  const size_t Ncells = mxGetNumberOfElements(mxVol);
  const size_t Ndofs = mxGetNumberOfElements(mxU0);
  const size_t Mspecies = Ndofs/Ncells;
  const size_t dsize = mxGetM(mxLData);  

  /* pointers to non-sparse objects */
  const double *vol = mxGetPr(mxVol);
  const double *ldata = mxGetPr(mxLData);
  const double *gdata = mxGetPr(mxGData);
  const double *sd_double = mxGetPr(mxSd);

  /* "typecast" from double to URDMEstate_t */
#ifdef UDS_
  /* make by the deterministic solver? */
  const URDMEstate_t *u0 = mxGetPr(mxU0);
  /* since URDMEstate_t == double */
#else
  const double *u0_double = mxGetPr(mxU0);
  URDMEstate_t *u0 = mxMalloc(Ndofs*sizeof(URDMEstate_t));
  for (int i = 0; i < Ndofs; i++) u0[i] = (URDMEstate_t)u0_double[i];
#endif

  /* typecast from double to int */
  int *sd = mxMalloc(Ncells*sizeof(int));
  for (int i = 0; i < Ncells; i++) sd[i] = (int)sd_double[i];

  /* parse solver arguments (K,I,S) */
  const double *K = mxGetPr(mxK);
  const size_t M1 = mxGetN(mxK);
  const double *I_double = mxGetPr(mxI);
  int *I = mxMalloc(3*M1*sizeof(int));
  for (int i = 0; i < 3*M1; i++) {
    I[i] = (int)I_double[i]-1;
    if (I[i] < 0 || Mspecies <= I[i])
      mexErrMsgTxt("Index out of bounds in inline propensity.");
  }
  int *prS;
  mwIndex *jcS;
  const bool emptyS = mxGetNumberOfElements(mxS) == 0;
  if (!emptyS) {
    jcS = mxGetJc(mxS);
    const double *prS_double = mxGetPr(mxS);
    const size_t nnzS = jcS[mxGetN(mxS)];
    /* typecast */
    prS = mxMalloc(nnzS*sizeof(int));
    for (int i = 0; i < nnzS; i++) prS[i] = (int)prS_double[i];
  }
  else
    /* empty S allowed */
    jcS = mxCalloc(M1+1,sizeof(mwIndex));

  /* result matrix R */
  plhs[0] = mxCreateDoubleMatrix(Mreactions,Ncells,mxREAL);
  double *R = mxGetPr(plhs[0]);

  /* fetch the propensities */
  PropensityFun *rfun;
  rfun = ALLOC_propensities(Mreactions-M1);

  for (int subvol = 0; subvol < Ncells; subvol++) {
    size_t j;
    for (j = 0; j < M1; j++)
      R[Mreactions*subvol+j] = 
	inlineProp(&u0[Mspecies*subvol],&K[j*3],&I[j*3],&prS[jcS[j]],
		   jcS[j+1]-jcS[j],vol[subvol],sd[subvol]);
    for (; j < Mreactions; j++)
      R[Mreactions*subvol+j] = 
	(*rfun[j-M1])(&u0[Mspecies*subvol],tt,vol[subvol],
		      &ldata[subvol*dsize],gdata,sd[subvol]);
  }

  /* deallocate */
  FREE_propensities(rfun);
  if (M1 > 0) {
    mxFree(I);
    if (!emptyS)
      mxFree(prS);
    else
      mxFree(jcS);
  }
  mxFree(sd);
#ifndef UDS_
  mxFree(u0);
#endif
}
/*----------------------------------------------------------------------*/
