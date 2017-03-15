/* mexnsm.c - Mex-interface NSM solver for use with URDME. */

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
  if (nrhs != 12 || nlhs != 1)
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
  const mxArray *mxSOLVEARGS = prhs[11];

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
  const double *vol = mxGetPr(mxVol);
  const double *ldata = mxGetPr(mxLData);
  const double *gdata = mxGetPr(mxGData);
  const double *sd_double = mxGetPr(mxSd);

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

  /* Parse solver arguments (K,I,S), if any */
  double *K;
  int *I = NULL,*prS;
  mwIndex *jcS;
  bool emptyS = true;
  size_t M1 = 0;
  if (!mxIsEmpty(mxSOLVEARGS)) {
    mxArray *mxK = NULL,*mxI = NULL,*mxS = NULL;
    const size_t nsolveargs = mxGetNumberOfElements(mxSOLVEARGS);
    if (nsolveargs&1 || !mxIsCell(mxSOLVEARGS))
      mexErrMsgTxt("Solver arguments must be an even length cell-vector.");
    for (int i = 0; i < nsolveargs; i += 2) {
      /* read property + value */
      mxArray *Str = mxGetCell(mxSOLVEARGS,i);
      char str[5];
      mxArray *Val = mxGetCell(mxSOLVEARGS,i+1);
      if (Str == NULL || Val == NULL)
	mexErrMsgTxt("Could not read property/value.");
      if (mxGetString(Str,str,2))
	mexErrMsgTxt("Could not parse property.");
      switch (str[0]) {
      case 'K':
	mxK = Val;
	M1 = mxGetN(mxK);
	break;
      case 'I':
	mxI = Val;
	break;
      case 'S':
	mxS = Val;
	emptyS = false;
	break;
      default:
	mexErrMsgTxt("Unrecognized property.");
      }
    }
    if (mxK == NULL || mxI == NULL)
      mexErrMsgTxt("One of (K,I) not parsed.");
    if (mxIsSparse(mxK) || !mxIsDouble(mxK) || mxGetM(mxK) != 3 || M1 > Mreactions ||
	mxIsSparse(mxI) || !mxIsDouble(mxI) || mxGetM(mxI) != 3 || mxGetN(mxI) != M1 || 
	(mxS != NULL && (!mxIsSparse(mxS) || !mxIsDouble(mxS) || mxGetN(mxS) != M1)))
      mexErrMsgTxt("Format mismatch in one of (K,I,S).");

    /* get pointer to K */
    K = mxGetPr(mxK);
    const double *I_double = mxGetPr(mxI);

    /* typecast I */
    I = mxMalloc(3*M1*sizeof(int));
    for (int i = 0; i < 3*M1; i++) {
      I[i] = (int)I_double[i]-1;
      if (I[i] < 0 || Mspecies <= I[i])
	mexErrMsgTxt("Index out of bounds in inline propensity.");
    }

    /* get sparse matrix S, if any */
    if (mxS != NULL) {
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
  }

  /* Typecast from double to int. */
  int *u0 = mxMalloc(Ndofs*sizeof(int));
  int *sd = mxMalloc(Ncells*sizeof(int));
  int *prN = mxMalloc(nnzN*sizeof(int));

  for (int i = 0; i < Ndofs; i++) u0[i] = (int)u0_double[i];
  for (int i = 0; i < Ncells; i++) sd[i] = (int)sd_double[i];
  for (int i = 0; i < nnzN; i++) prN[i] = (int)prN_double[i];

  /* Set up result matrix U. */
  int *U = mxMalloc(Ndofs*tlen*sizeof(int));

  /* report */
  const int report_level = (int)*mxGetPr(mxREPORT);

  /* seed */
  if (!mxIsNaN(*mxGetPr(mxSEED))) srand48((long)*mxGetPr(mxSEED));

  /* fetch the propensities, execute the Master call to nsm(), and
     next get rid of the propensities */
  PropensityFun *rfun;
  rfun = ALLOC_propensities(Mreactions-M1);
  nsm(rfun,
      (const size_t *)irD,(const size_t *)jcD,prD,
      u0,
      (const size_t *)irN,(const size_t *)jcN,prN,
      (const size_t *)irG,(const size_t  *)jcG,
      tspan,tlen,
      U,vol,ldata,gdata,sd,
      Ncells,Mspecies,Mreactions,
      dsize,
      report_level,
      K,I,(const size_t *)jcS,prS,M1
      );
  FREE_propensities(rfun);

  /* Deallocate. */
  if (M1 > 0) {
    mxFree(I);
    if (!emptyS)
      mxFree(prS);
    else
      mxFree(jcS);
  }
  mxFree(prN);
  mxFree(sd);
  mxFree(u0);

  /* Put result in plhs[0] and typecast from int to double. */
  plhs[0] = mxCreateDoubleMatrix(Ndofs,tlen,mxREAL);
  double *U_out = mxGetPr(plhs[0]);
  for (int i = 0; i < Ndofs*tlen; i++) U_out[i] = (double)U[i];
  mxFree(U);
}
/*----------------------------------------------------------------------*/
