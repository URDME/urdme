/* mexrhs.c - Mex-interface of propensity evaluation for use with the
   UDS solver. */

/* S. Engblom 2024-05-13 (ldata_time, gdata_time) */
/* S. Engblom 2019-11-27 (Revision, inline propensities) */
/* S. Engblom 2019-11-06 (Revision, now using URDMEstate_t) */
/* S. Engblom 2017-02-24 */

#include <stdlib.h>

#include "mex.h"
#include "matrix.h"

#include "inline.h"
#include "propensities.h"
#include "report.h"

/*----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /* check syntax */
  if (nrhs != 13) mexErrMsgTxt("Wrong number of arguments.");
  if (*(unsigned *)mxGetData(prhs[0]) != MEXHASH)
    PERROR("Mex hashkey mismatch.");

  /* load arguments */
  const mxArray *mxTspan = prhs[1];
  const mxArray *mxU0 = prhs[2];
  const size_t Mreactions = (size_t)mxGetScalar(prhs[3]);
  const mxArray *mxVol = prhs[4];
  const mxArray *mxLData = prhs[5];
  const mxArray *mxGData = prhs[6];
  const mxArray *mxLDataTime = prhs[7];
  const mxArray *mxGDataTime = prhs[8];
  const mxArray *mxSd = prhs[9];
  const mxArray *mxK = prhs[10];
  const mxArray *mxI = prhs[11];
  const mxArray *mxS = prhs[12];

  /* get problem dimensions */
  const size_t Nreplicas = mxGetNumberOfDimensions(mxU0) == 3
    ? mxGetDimensions(mxU0)[2] : 1;
  const size_t tlen = mxGetNumberOfElements(mxTspan);
  const size_t Ncells = mxGetNumberOfElements(mxVol);
  const size_t Ndofs = mxGetNumberOfElements(mxU0)/tlen/Nreplicas;
  const size_t Mspecies = Ndofs/Ncells;
  const size_t ldsize = mxGetM(mxLData);
  const size_t ldtsize = mxGetM(mxLDataTime);
  
  /* pointers to non-sparse objects */
  const double *tspan = mxGetPr(mxTspan);
  const double *vol = mxGetPr(mxVol);
  const double *ldata = mxGetPr(mxLData);
  const double *gdata = mxGetPr(mxGData);
  const double *ldata_time = mxGetPr(mxLDataTime);
  const double *gdata_time = mxGetPr(mxGDataTime);
  const double *sd_double = mxGetPr(mxSd);

  /* "typecast" from double to URDMEstate_t */
#ifdef UDS_
  /* make by the deterministic solver? */
  const URDMEstate_t *u0 = mxGetPr(mxU0);
  /* since URDMEstate_t == double */
#else
  const double *u0_double = mxGetPr(mxU0);
  URDMEstate_t *u0 = mxMalloc(mxGetNumberOfElements(mxU0)*sizeof(URDMEstate_t));
  for (int i = 0; i < mxGetNumberOfElements(mxU0); i++)
    u0[i] = (URDMEstate_t)u0_double[i];
#endif

  /* typecast from double to int */
  int *sd = mxMalloc(Ncells*sizeof(int));
  for (int i = 0; i < Ncells; i++) sd[i] = (int)sd_double[i];

  /* Parse inline propensities (K,I,S) - assume empty to begin with */
  double *K = NULL;
  int *I = NULL,*prS = NULL;
  mwIndex *jcS = NULL;
  const bool emptyS = mxIsEmpty(mxS);
  size_t M1 = 0;

  /* both of {K,I} empty is a common case */
  if (!mxIsEmpty(mxK) || !mxIsEmpty(mxI)) {
    /* get pointer to K */
    K = mxGetPr(mxK);
    M1 = mxGetN(mxK);

    /* typecast I */
    const double *I_double = mxGetPr(mxI);
    I = mxMalloc(3*M1*sizeof(int));
    for (int i = 0; i < 3*M1; i++)
      I[i] = (int)I_double[i]-1;
  }

  /* get sparse matrix S, if any */
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
  if (Nreplicas == 1)
    plhs[0] = mxCreateDoubleMatrix(Mreactions*Ncells,tlen,mxREAL);
  else {
    const mwSize dims[] = {Mreactions*Ncells,tlen,Nreplicas};
    plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
  }
  double *R = mxGetPr(plhs[0]);

  /* fetch the propensities */
  PropensityFun *rfun;
  rfun = ALLOC_propensities(Mreactions-M1);

  /* immediate loops */
  for (int k = 0; k < Nreplicas; k++)
    for (int j = 0; j < tlen; j++)
      for (int subvol = 0; subvol < Ncells; subvol++) {
	size_t i;
	for (i = 0; i < M1; i++)
	  R[tlen*Ncells*Mreactions*k+Ncells*Mreactions*j+Mreactions*subvol+i] = 
	    inlineProp(&u0[tlen*Ncells*Mspecies*k+
			   Ncells*Mspecies*j+
			   Mspecies*subvol],
		       &K[i*3],&I[i*3],&prS[jcS[i]],
		       jcS[i+1]-jcS[i],vol[subvol],sd[subvol]);
	for (; i < Mreactions; i++)
	  R[tlen*Ncells*Mreactions*k+Ncells*Mreactions*j+Mreactions*subvol+i] = 
	    (*rfun[i-M1])(&u0[tlen*Ncells*Mspecies*k+
			      Ncells*Mspecies*j+
			      Mspecies*subvol],
			  tspan[j],vol[subvol],
			  &ldata[subvol*ldsize],gdata,
			  &ldata_time[subvol*ldtsize],gdata_time,
			  sd[subvol]);
      }
  FREE_propensities(rfun);

  /* deallocate */
  if (M1 > 0) {
    mxFree(I);
    mxFree(prS);
    if (emptyS) mxFree(jcS);
  }
  mxFree(sd);
#ifndef UDS_
  mxFree(u0);
#endif
}
/*----------------------------------------------------------------------*/
