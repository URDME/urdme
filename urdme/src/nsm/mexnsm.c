/* mexnsm.c - Mex-interface NSM solver for use with URDME. */

/* S. Engblom 2024-05-08 (data_time, ldata_time, gdata_time) */
/* S. Engblom 2019-11-27 (Revision, inline propensities) */
/* S. Engblom 2019-11-12 (multiple seeds) */
/* S. Engblom 2018-02-10 (Nreplicas syntax) */
/* S. Engblom 2017-02-16 (Major revision, URDME 1.3, Comsol 5) */
/* J. Cullhed 2008-08-04 (mexrdme.c) . */

#include <stdlib.h>

#include "mex.h"
#include "matrix.h"

#include "propensities.h"
#include "nsm.h"
#include "report.h"

/*----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /* check syntax */
  if (nrhs != 19 || nlhs != 1)
    mexErrMsgTxt("Wrong number of arguments.");
  if (*(unsigned *)mxGetData(prhs[0]) != MEXHASH)
    PERROR("Mex hashkey mismatch.");

  /* load arguments */
  const mxArray *mxTspan = prhs[1];
  const mxArray *mxU0 = prhs[2];
  const mxArray *mxD = prhs[3];
  const mxArray *mxN = prhs[4];
  const mxArray *mxG = prhs[5];
  const mxArray *mxVol = prhs[6];
  const mxArray *mxLData = prhs[7];
  const mxArray *mxGData = prhs[8];
  const mxArray *mxDataTime = prhs[9];
  const mxArray *mxLDataTime = prhs[10];
  const mxArray *mxGDataTime = prhs[11];
  const mxArray *mxSd = prhs[12];
  const mxArray *mxREPORT = prhs[13];
  const mxArray *mxSEED = prhs[14];
  const mxArray *mxK = prhs[15];
  const mxArray *mxI = prhs[16];
  const mxArray *mxS = prhs[17];
  const mxArray *mxSOLVEARGS = prhs[18];

  /* get dimensions */
  const size_t tlen = mxGetNumberOfElements(mxTspan);
  const size_t Ncells = mxGetNumberOfElements(mxVol);
  const size_t Mspecies = mxGetM(mxN);
  const size_t Mreactions = mxGetN(mxN);
  const size_t ldsize = mxGetM(mxLData);
  const size_t ldtsize = mxGetM(mxLDataTime);
  const size_t gdtsize = mxGetM(mxGDataTime);
  const size_t dtlen =  mxGetN(mxGDataTime);
  const size_t Ndofs = Ncells*Mspecies;

  /* Get pointers to non-sparse objects. */
  const double *tspan = mxGetPr(mxTspan);
  const double *u0_double = mxGetPr(mxU0);
  const size_t Nreplicas = mxGetNumberOfDimensions(mxU0) == 3
    ? mxGetDimensions(mxU0)[2] : 1;
  const double *vol = mxGetPr(mxVol);
  const double *ldata = mxGetPr(mxLData);
  const double *gdata = mxGetPr(mxGData);
  const double *data_time = mxGetPr(mxDataTime);
  const double *ldata_time = mxGetPr(mxLDataTime);
  const double *gdata_time = mxGetPr(mxGDataTime);
  const double *sd_double = mxGetPr(mxSd);
  const double *seed_double = mxGetPr(mxSEED);

  /* Get sparse matrix D. */
  const mwIndex *jcD = (const mwIndex *)mxGetJc(mxD);
  const mwIndex *irD = (const mwIndex *)mxGetIr(mxD);
  const double *prD = mxGetPr(mxD);

  /* Get sparse matrix N. */
  const mwIndex *jcN = (const mwIndex *)mxGetJc(mxN);
  const mwIndex *irN = (const mwIndex *)mxGetIr(mxN);
  const double *prN_double = mxGetPr(mxN);
  const mwIndex nnzN = jcN[mxGetN(mxN)];

  /* Get sparse matrix G. */
  const mwIndex *jcG = (const mwIndex *)mxGetJc(mxG);
  const mwIndex *irG = (const mwIndex *)mxGetIr(mxG);

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

  /* Parse solver arguments */
  if (!mxIsEmpty(mxSOLVEARGS))
    mexErrMsgTxt("NSM does not accept any solver arguments");
 
  /* Typecast from double to int. */
  int *u0 = mxMalloc(Ndofs*Nreplicas*sizeof(int));
  int *sd = mxMalloc(Ncells*sizeof(int));
  int *prN = mxMalloc(nnzN*sizeof(int));

  for (int i = 0; i < Ndofs*Nreplicas; i++) u0[i] = (int)u0_double[i];
  for (int i = 0; i < Ncells; i++) sd[i] = (int)sd_double[i];
  for (int i = 0; i < nnzN; i++) prN[i] = (int)prN_double[i];

  /* Set up result matrix U. */
  int *U = mxMalloc(Ndofs*tlen*Nreplicas*sizeof(int));

  /* report */
  const int report_level = (int)*mxGetPr(mxREPORT);

  /* seed */
  long *seed_long = mxMalloc(Nreplicas*sizeof(long));
  for (int i = 0,j = 0,
	 jinc = mxGetNumberOfElements(mxSEED) == Nreplicas;
       i < Nreplicas; i++, j += jinc)
    seed_long[i] = (long)seed_double[j];

  /* fetch the propensities, execute the Master call to nsm(), and
     next get rid of the propensities */
  PropensityFun *rfun;
  rfun = ALLOC_propensities(Mreactions-M1);
  nsm(rfun,
      (const size_t *)irD,(const size_t *)jcD,prD,
      u0,
      (const size_t *)irN,(const size_t *)jcN,prN,
      (const size_t *)irG,(const size_t  *)jcG,
      tspan,tlen,Nreplicas,
      U,vol,ldata,gdata,data_time,ldata_time,gdata_time,sd,
      Ncells,Mspecies,Mreactions,
      ldsize,ldtsize,gdtsize,dtlen,
      report_level,seed_long,
      K,I,(const size_t *)jcS,prS,M1
      );
  FREE_propensities(rfun);

  /* Deallocate. */
  if (M1 > 0) {
    mxFree(I);
    mxFree(prS);
    if (emptyS) mxFree(jcS);
  }
  mxFree(seed_long);
  mxFree(prN);
  mxFree(sd);
  mxFree(u0);

  /* Put result in plhs[0] and typecast from int to double. */
  if (Nreplicas == 1)
    plhs[0] = mxCreateDoubleMatrix(Ndofs,tlen,mxREAL);
  else {
    const mwSize dims[] = {Ndofs,tlen,Nreplicas};
    plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
  }
  double *U_out = mxGetPr(plhs[0]);
  for (int i = 0; i < Ndofs*tlen*Nreplicas; i++) U_out[i] = (double)U[i];
  mxFree(U);
}
/*----------------------------------------------------------------------*/
