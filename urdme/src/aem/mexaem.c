/* mexaem.c - Mex-interface AEM solver for use with URDME. */

/* S. Engblom 2019-11-12 (Nreplicas syntax) */
/* S. Engblom 2017-02-16 (Major revision, URDME 1.3, Comsol 5) */
/* J. Cullhed 2008-08-04 (mexrdme.c) . */

#include <stdlib.h>

#include "mex.h"
#include "matrix.h"

#include "propensities.h"
#include "aem.h"
#include "report.h"

/*----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /* check syntax */
  if (nrhs != 15 || nlhs != 1)
    mexErrMsgTxt("Wrong number of arguments.");

  /* load arguments */
  const mxArray *mxTspan = prhs[0];
  const mxArray *mxU0 = prhs[1];
  const mxArray *mxD = prhs[2];
  const mxArray *mxN = prhs[3];
  const mxArray *mxG = prhs[4];
  const mxArray *mxVol = prhs[5];
  const mxArray *mxLData = prhs[6];
  const mxArray *mxGData = prhs[7];
  const mxArray *mxSd = prhs[8];
  const mxArray *mxREPORT = prhs[9];
  const mxArray *mxSEED = prhs[10];
  const mxArray *mxK = prhs[11];
  const mxArray *mxI = prhs[12];
  const mxArray *mxS = prhs[13];
  const mxArray *mxSOLVEARGS = prhs[14];

  /* get dimensions */
  const size_t tlen = mxGetNumberOfElements(mxTspan);
  const size_t Ncells = mxGetNumberOfElements(mxVol);
  const size_t Mspecies = mxGetM(mxN);
  const size_t Mreactions = mxGetN(mxN);
  const size_t dsize = mxGetM(mxLData);
  const size_t Ndofs = Ncells*Mspecies;

  /* Get pointers to non-sparse objects. */
  const double *tspan = mxGetPr(mxTspan);
  const double *u0_double = mxGetPr(mxU0);
  const size_t Nreplicas = mxGetNumberOfDimensions(mxU0) == 3
    ? mxGetDimensions(mxU0)[2] : 1;
  const double *vol = mxGetPr(mxVol);
  const double *ldata = mxGetPr(mxLData);
  const double *gdata = mxGetPr(mxGData);
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
    mexErrMsgTxt("AEM does not accept any solver arguments");

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
  unsigned *seed_uint = mxMalloc(Nreplicas*sizeof(unsigned));
  for (int i = 0,j = 0,
	 jinc = mxGetNumberOfElements(mxSEED) == Nreplicas;
       i < Nreplicas; i++, j += jinc)
    seed_uint[i] = (unsigned)seed_double[j];
  
  /* fetch the propensities, execute the Master call to aem(), and
     next get rid of the propensities */
  PropensityFun *rfun;
  rfun = ALLOC_propensities(Mreactions-M1);
  aem(rfun,
      (const size_t *)irD,(const size_t *)jcD,prD,
      u0,
      (const size_t *)irN,(const size_t *)jcN,prN,
      (const size_t *)irG,(const size_t  *)jcG,
      tspan,tlen,Nreplicas,
      U,vol,ldata,gdata,sd,
      Ncells,Mspecies,Mreactions,
      dsize,
      report_level,seed_uint,
      K,I,(const size_t *)jcS,prS,M1);
  FREE_propensities(rfun);

  /* Deallocate. */
  if (M1 > 0) {
    mxFree(I);
    mxFree(prS);
    if (emptyS) mxFree(jcS);
  }
  mxFree(seed_uint);
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
