/* mexrhs.c - Mex-interface of Jacobian of propensities. */

/* S. Engblom 2020-02-21 (based on mexrhs.c) */

#include <stdlib.h>
#include <string.h>

#include "mex.h"
#include "matrix.h"

#include "inline.h"
#include "propensities.h"
#include "report.h"

/*----------------------------------------------------------------------*/
double Diff_inlineProp(const URDMEstate_t *xx,const double *k,const int *i,
		       const int *prS,size_t nS,int species,double vol,int sd)
/* Derivative of inline propensity wrt species. */
{
  double res = 0.0;

  /* propensity turned off? */
  if (sd == 0) return 0.0;
  for (size_t j = 0; j < nS; j++) if (sd == prS[j]) return 0.0;

  /* main formula, check for 0 rate constants to avoid index out of
     bounds */
  if (k[0] != 0.0) {
    if (i[0] != i[1]) {
      if (i[0] == species)
	res += k[0]*xx[i[1]]/vol;
      else
	res += k[0]*xx[i[0]]/vol;
    }
    else if (i[0] == species)
      res += k[0]*(2.0*xx[i[0]]-1.0)/(2.0*vol);
  }
  if (k[1] != 0.0 && i[2] == species) res += k[1];

  return res;
}
/*----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /* check syntax */
  if (nrhs != 14) mexErrMsgTxt("Wrong number of arguments.");
  if (*(unsigned *)mxGetData(prhs[0]) != MEXHASH)
    PERROR("Mex hashkey mismatch.");

  /* load arguments */
  const double tt = mxGetScalar(prhs[1]);
  const mxArray *mxU0 = prhs[2];
  const size_t Mreactions = (size_t)mxGetScalar(prhs[3]);
  const mxArray *mxG = prhs[4];
  const mxArray *mxVol = prhs[5];
  const mxArray *mxLData = prhs[6];
  const mxArray *mxGData = prhs[7];
  const mxArray *mxLDataTime = prhs[8];
  const mxArray *mxGDataTime = prhs[9];
  const mxArray *mxSd = prhs[10];
  const mxArray *mxK = prhs[11];
  const mxArray *mxI = prhs[12];
  const mxArray *mxS = prhs[13];

  /* get problem dimensions */
  const size_t Ncells = mxGetNumberOfElements(mxVol);
  const size_t Ndofs = mxGetNumberOfElements(mxU0);
  const size_t Mspecies = Ndofs/Ncells;
  const size_t ldsize = mxGetM(mxLData);
  const size_t ldtsize = mxGetM(mxLDataTime);
  
  /* pointers to non-sparse objects */
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
  URDMEstate_t *u0 = mxMalloc(Ndofs*sizeof(URDMEstate_t));
  for (int i = 0; i < Ndofs; i++) u0[i] = (URDMEstate_t)u0_double[i];
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

  /* dependency matrix G - only the diffusion part is used since this
     tells us the dependency of reactions upon species */
  const mwIndex *jcG = mxGetJc(mxG);
  const mwIndex *irG = mxGetIr(mxG);
  const size_t nnzG = jcG[Mspecies]-jcG[0];

  /* working memory to hold Jacobian of compiled propensities */
  double *val = NULL;
#ifdef JAC_
  if (M1 < Mreactions)
    val = mxMalloc(nnzG*sizeof(double));
#else
    if (M1 < Mreactions)
      mexErrMsgTxt("Attempt to compute Jacobian of propensities with "	\
		   "no symbolic differentiation available. Compile " \
		   "with JAC_ #defined and ensure that the propensity " \
		   "file defines GET_jacobian().");
#endif

  /* allocate Jacobian sparse matrix JAC */
  plhs[0] = mxCreateSparse(Mreactions*Ncells,Mspecies*Ncells,
			   nnzG*Ncells,mxREAL);
  mwIndex *jcJAC = mxGetJc(plhs[0]);
  mwIndex *irJAC = mxGetIr(plhs[0]);
  double *prJAC = mxGetPr(plhs[0]);

  /* sparsity pattern deduced from G */
  memcpy(&jcJAC[0],&jcG[0],(Mspecies+1)*sizeof(mwIndex));
  memcpy(&irJAC[0],&irG[jcG[0]],nnzG*sizeof(mwIndex));
  /* the global sparsity pattern is a simple block-diagonal shift of
     the pattern in the first subvolume: so loop over subvolumes and
     repeat pattern from previous block while adding an offset */
  for (int subvol = 1; subvol < Ncells; subvol++) {
    /* columns */
    memcpy(&jcJAC[subvol*Mspecies+1],&jcJAC[(subvol-1)*Mspecies+1],
	   Mspecies*sizeof(mwIndex));
    for (int j = 0; j < Mspecies; j++)
      jcJAC[j+subvol*Mspecies+1] += nnzG;
    /* rows */
    memcpy(&irJAC[subvol*nnzG],&irJAC[(subvol-1)*nnzG],
	   nnzG*sizeof(mwIndex));
    for (int i = 0; i < nnzG; i++)
      irJAC[i+subvol*nnzG] += Mreactions;
  }

  /* once the parsity pattern is in place, perform the
     differentiation */
  for (int subvol = 0; subvol < Ncells; subvol++) {
    /* get the Jacobian from all compiled propensities */
    int valix = 0;
#ifdef JAC_
    GET_jacobian(val,&u0[Mspecies*subvol],tt,vol[subvol],
		 &ldata[subvol*ldsize],gdata,
		 &ldata_time[subvol*ldtsize],gdata_time,
		 sd[subvol]);
#endif

    /* loop over columns/species in this subvolume... */
    for (int j = 0; j < Mspecies; j++)
      /* ...and loop over nonzero rows/reactions */
      for (int i = jcJAC[j+subvol*Mspecies]; i < jcJAC[j+subvol*Mspecies+1];
	   i++) {
	/* after figuring out the reaction, compute the derivative */
	const int reaction = irJAC[i]%Mreactions;
	if (reaction < M1)
	  prJAC[i] =
	    Diff_inlineProp(&u0[Mspecies*subvol],
			    &K[reaction*3],&I[reaction*3],&prS[jcS[reaction]],
			    jcS[reaction+1]-jcS[reaction],j,
			    vol[subvol],sd[subvol]);
	else
	  prJAC[i] = val[valix++];
      }
  }

  /* deallocate */
  if (M1 < Mreactions)
    mxFree(val);
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
