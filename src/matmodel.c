/* A. Hellander and B. Drawert. */
#include <string.h>
#include "matmodel.h"

// This is necessary as 'parameters' is extern in "propensities.h" (where it is declared).
//      It must be defined (set to a value) in one and only one '.o' file
#include "propensities.h"
double *parameters = NULL;

/*
 
 URDME_LIBMAT is defined at compile time, 
 and dictates wether URMDE will use Mathworks' MATLAB libraries
 for handling .mat files, or our own implmentation.
 
 The Matlab libraries are currently needed for 
 sparse output.  
 
 */

#ifdef URDME_LIBMAT
#include "read_matfile.h"
#else
#include "mat.h" 
#include "mex.h"
#include "matrix.h"
#endif

/* 
 
   Read model input file (.mat format) and initialize urdme model struct. 
   If any of the required fields (for the nsm solver) is missing, 
   it returns a NULL-pointer.
   
*/



urdme_model *read_model(char *file)
{
	
	int i,Ndofs;
	
	mxArray *D,*N,*G,*vol,*u0,*tspan,*data,*sd;
	
	urdme_model *model;
	model = (urdme_model *)malloc(sizeof(urdme_model));
	
	/* Open mat-file (read-only) */
	MATFile *input_file;
	input_file = matOpen(file,"r"); 
	
	if (input_file == NULL){
		printf("Failed to open mat-file.\n");
		return NULL;	
	}
	
	/* Get D-matrix. */
	D = matGetVariable(input_file, "D");
	if (D == NULL){
		printf("The diffusion matrix D is missing in the model file.\n");
		return NULL;
	}
	
    Ndofs  = mxGetN(D); 
    size_t *mxirD;
    size_t *mxjcD;
    double *mxprD;
    mxirD = mxGetIr(D);
    mxjcD = mxGetJc(D);
    mxprD = mxGetPr(D);
    int nnzD = mxGetNzmax(D);
    
    model->jcD = (size_t *)malloc((Ndofs+1)*sizeof(size_t));
    memcpy(model->jcD,mxjcD,(Ndofs+1)*sizeof(size_t));
    model->irD = (size_t*)malloc(nnzD*sizeof(size_t));
    memcpy(model->irD,mxirD,nnzD*sizeof(size_t));
    model->prD = (double *)malloc(nnzD*sizeof(double));
    memcpy(model->prD,mxprD,nnzD*sizeof(double));
    mxDestroyArray(D);
    
	/* Stoichiometry matrix */
	N = matGetVariable(input_file, "N");
    if (N == NULL){
		printf("The stoichiometry matrix is missing in the model file.\n");
		return NULL;
	}
  
	/* Typecast to int */
	double *tempN;
	tempN = mxGetPr(N);
	int nnzN = (int) mxGetNzmax(N);
	int *prN;
	prN = (int*)malloc(nnzN*sizeof(int));
	for (i=0;i<nnzN;i++)
		prN[i] = (int) tempN[i];
	
    
    model->irN = (size_t *)malloc(nnzN*sizeof(size_t));
    memcpy(model->irN,(size_t *)mxGetIr(N),nnzN*sizeof(size_t));
    model->jcN = (size_t *)malloc((mxGetN(N)+1)*sizeof(size_t));
    memcpy(model->jcN,(size_t *)mxGetJc(N),(mxGetN(N)+1)*sizeof(size_t));
    
	model->prN = prN;
	model->Mspecies   = (int) mxGetM(N);
	model->Mreactions = (int) mxGetN(N);
	model->Ncells=Ndofs/model->Mspecies;
    
    mxDestroyArray(N);
	
    /* Volume vector */
	vol = matGetVariable(input_file,"vol");
	if (vol == NULL){
		printf("The volume vector is missing in the model file.\n");
		return NULL;
	}
    model->vol = (double *)malloc(model->Ncells*sizeof(double));
    memcpy(model->vol,(double *)mxGetPr(vol),model->Ncells*sizeof(double));
    mxDestroyArray(vol);
    
	/* Dependency graph */
	G = matGetVariable(input_file, "G");
    if (G == NULL){
		printf("The dependency graph (G) is missing in the model file.\n");
		return NULL;
	}
    int nnzG = mxGetNzmax(G);
    model->irG = (size_t *)malloc(nnzG*sizeof(size_t));
    memcpy(model->irG,(size_t *)mxGetIr(G),nnzG*sizeof(size_t));
	model->jcG = (size_t *)malloc((mxGetN(G)+1)*sizeof(size_t));
    memcpy(model->jcG,(size_t *)mxGetJc(G),(mxGetN(G)+1)*sizeof(size_t));
	mxDestroyArray(G);
	
	/* initial condition */
	u0 = matGetVariable(input_file, "u0");
	if (u0 == NULL){
		printf("Initial condition (u0) is missing in the model file.\n");
		return NULL;
	}
	/* Typecast */
	
	int *u0int;
	u0int = (int*)malloc(Ndofs*sizeof(int));
	double *u0temp;
	u0temp = mxGetPr(u0);
	for (i=0;i<Ndofs;i++)
		u0int[i] = (int) u0temp[i];
	model->u0 = u0int;
    mxDestroyArray(u0);

    /* time span */
	tspan = matGetVariable(input_file, "tspan");
	if (tspan == NULL){
		printf("Time span (tspan) is missing in the model file.\n");
		return NULL;
	}
	model->tlen = mxGetNumberOfElements(tspan);
    model->tspan = (double *)malloc(model->tlen*sizeof(double));
    memcpy(model->tspan,(double *)mxGetPr(tspan),model->tlen*sizeof(double));
    mxDestroyArray(tspan);
    
	/* Subdomain vector */
	sd = matGetVariable(input_file, "sd");
	if (sd == NULL){
		printf("Subdomain vector (sd) is missing in the model file.\n");
		return NULL;
	}
	
	/* typecast */
	int *sdint = (int*)malloc(model->Ncells*sizeof(int));
	double *sdtemp;
	sdtemp = mxGetPr(sd);
	for (i=0;i<model->Ncells;i++)
		sdint[i]=(int) sdtemp[i];
    
	model->sd = sdint;
    mxDestroyArray(sd);

    
	/* data matrix */
	data = matGetVariable(input_file, "data");
	if (data == NULL){
		printf("Data matrix (data) is missing in the model file.\n");
		return NULL;
	}
	
	model->dsize = mxGetM(data);
    model->data = (double *)malloc(model->dsize*model->Ncells*sizeof(double));
    memcpy(model->data,(double *)mxGetPr(data),model->dsize*model->Ncells*sizeof(double));
	//model->data  = mxGetPr(data);
    mxDestroyArray(data);
	
    /* Maximum number of solutions defaults to one. */
    model->nsolmax = 1;
    
	matClose(input_file);	
	
	return model;
	
	
}

//--------------------------------------------
void read_solution(urdme_model *model, char*file){
	mxArray *tspan,*U;
	MATFile *input_file;
	/* Open mat-file (read-only) */
	input_file = matOpen(file,"r"); 
	if (input_file == NULL){
		printf("Failed to open mat-file '%s'.\n",file);
		return;	
	}
    printf("read in '%s'\n'",file);
	
	/* Get U-matrix. */
    init_sol(model,1);
    U = matGetVariable(input_file, "U");
	if (U == NULL){
        printf("No U solution variable in mat-file '%s'\n.",file);
        return;
    }
    int Usz = mxGetNumberOfElements(U);
	model->U[0] =(int *)malloc(Usz*sizeof(int));
    double*Utmp = mxGetPr(U);
	model->nsol = 1;
    int i;
    for(i=0;i<Usz;i++){
        model->U[0][i] = (int) Utmp[i];
    }

	/* time span (optional) */
	tspan = matGetVariable(input_file, "tspan");
	if (tspan != NULL){
        model->tspan = mxGetPr(tspan);
        model->tlen = mxGetNumberOfElements(tspan);
	}
}


/* Utility function. Initialize the solution field of the urdme_model struct. */
void init_sol(urdme_model *model, int nsolmax)
{
	int i;
	model->nsolmax = nsolmax;
	model->U = (int**)malloc(model->nsolmax*sizeof(int *));
	for (i=0;i<nsolmax;i++)
		model->U[i]=NULL;
	model->nsol = 0;
	
}

/* Free model. Always call this destructor before exiting. */
int destroy_model(urdme_model *model)
{
    
    int i;
    
    /* D-matrix */
    free(model->irD);
    free(model->jcD);
    free(model->prD);
    
    /* Volume vector */
    free(model->vol);
    
    /* Dependency graph */
    free(model->jcG);
    free(model->irG);
    
    /* Stoichiometry matrix */
    free(model->irN);
    free(model->jcN);
    free(model->prN);
	
    free(model->tspan);
    free(model->u0);
	free(model->sd);
    free(model->data);
	free(model->U);
    
    for (i=0; i<model->num_extra_args; i++)
      if (model->extra_args[i]!=NULL)
        free(model->extra_args[i]);
    
	free(model);
    
	return 0;
}

/* 
   Print trajectory attached to urdme_model to a .mat file.
   
   The output trajectory will be stored in a variable called "U", that can
   subsequently be loaded into the Matlab workspace for postprocessing. 
 
 */

int dump_results(urdme_model* model, char *filename, char *type){
	
	int Ndofs,tlen,nsol,i,j;
	int *U;

	MATFile *output;
	
	Ndofs = model->Ncells*model->Mspecies;
	tlen  = model->tlen;
	nsol  = model->nsol;

    mxArray *Uout; 
	U = model->U[0];

	
#ifndef URDME_OUTPUT_SPARSE
	Uout = mxCreateDoubleMatrix(Ndofs,tlen,mxREAL);
	double *data;
	data = mxGetPr(Uout);
	/* Dense output data */
	for (i=0;i<Ndofs*model->tlen;i++){
		data[i]=(double) U[i];
	}
    /* 
       
       If sparse output, we write the matrix on dictionary of keys (DOK) format.
       The application using this matrix will then have to convert to whatever format needed. 
       the keys are (iU,jU) and values sU. To instantiate a sparse matrix in Matlab, load the output file and do
       
       >> U = sparse(iU,jU,sU,mU,nU);
     
     */
#else	
	/* Count number of non-zeros */
	int nnz = 0,nnz_col;
	for (i=0; i<Ndofs*tlen; i++) {
		if (U[i]>0.0)
			nnz++;
	}
	
	/*Uout = mxCreateSparse(Ndofs,tlen,nnz,mxREAL);
	if (Uout == NULL) {
		printf("Fatal error. Failed to create output matrix.");
		exit(-1);
	}
	
	size_t *jc = mxGetJc(Uout);
	size_t *ir = mxGetIr(Uout);
	double *pr = mxGetPr(Uout);*/
    
    mxArray* iU = mxCreateDoubleMatrix(nnz,1,mxREAL);
    mxArray* jU = mxCreateDoubleMatrix(nnz,1,mxREAL);
    mxArray* sU = mxCreateDoubleMatrix(nnz,1,mxREAL);
    
    // Dimesions of the matrix, mU (row) and nU(col). In house matlib does not have mxCreateScalar.
    mxArray* mU = mxCreateDoubleMatrix(1,1,mxREAL);
    double *dim;
    dim = mxGetPr(mU);
    dim[0] = (double)Ndofs;
    mxArray* nU = mxCreateDoubleMatrix(1,1,mxREAL);
    dim = mxGetPr(nU);
    dim[0] = (double)tlen;
    
    double *iUdata = mxGetPr(iU);
    double *jUdata =  mxGetPr(jU);
    double *sUdata = mxGetPr(sU);
    
    nnz = 0;
    for (j=0; j<tlen; j++){
        for (i=0;i<Ndofs;i++){
            if (U[j*Ndofs+i]>0.0){
                 /* NOTE THE +1 HERE, MATLAB SPECIFIC INDEXING. */
                iUdata[nnz] = i+1;
                jUdata[nnz] = j+1;
                sUdata[nnz] = (double)U[j*Ndofs+i];
                nnz++;
            }
        }
    }
	
	/*nnz  = 0; jc[0]=0;
	for (j=0; j<tlen; j++) {
		for (i=0;i<Ndofs;i++) {
			
			if (U[j*Ndofs+i]>0) {
				pr[nnz] = (double)U[j*Ndofs+i];
				ir[nnz] = i;
				nnz++;
				nnz_col++;
			}
			
		}
		jc[j+1]=nnz;
	}*/
#endif
	
	
    mxArray *Tspanout; 
	Tspanout = mxCreateDoubleMatrix(1, model->tlen,mxREAL);
	double *tdata;
	tdata = mxGetPr(Tspanout);
		
	for(i=0;i<model->tlen;i++)
		tdata[i]=model->tspan[i];
	
	
#ifdef URDME_LIBMAT 
	output = matOpen(filename,"w");  // using built in URDME mat read/write
#else
	output = matOpen(filename,"w7.3"); // Possibly large files 
#endif

    if (output == NULL){
		printf("Error: Could not write to output file: %s\n",filename);
		return(-1);
	}
	
#ifndef URDME_OUTPUT_SPARSE
	matPutVariable(output,"U",Uout);
	matPutVariable(output,"tspan",Tspanout);
#else
    matPutVariable(output,"iU",iU);
    matPutVariable(output,"jU",jU);
    matPutVariable(output,"sU",sU);
    matPutVariable(output,"mU",mU);
    matPutVariable(output,"nU",nU);
	matPutVariable(output,"tspan",Tspanout);
#endif
    
	matClose(output);
    
#ifndef URDME_OUTPUT_SPARSE
	mxDestroyArray(Uout);
    mxDestroyArray(Tspanout);
#else
    mxDestroyArray(iU);
    mxDestroyArray(jU);
    mxDestroyArray(sU);
    mxDestroyArray(Tspanout);
    mxDestroyArray(mU);
    mxDestroyArray(nU);
#endif
    
	
	return 1;
}

/*
    Print all data from model to stdout, for debugging
*/
void debug_print_model(urdme_model* model){
    int i,j;
	/* Constants */
	printf("Mspecies=%i\n",model->Mspecies);
	printf("Mreactions=%i\n",model->Mreactions);
	printf("Ncells=%i\n",model->Ncells);

    int Ndofs = model->Ncells * model->Mspecies;
    
    for (i=0; i<Ndofs; i++) {
        for (j=model->jcD[i]; j<model->jcD[i+1]; j++){
            printf("D: %i\t%i\t%f\n",j,(int)model->irD[j],model->prD[j]);
        }
    }
    printf("end D\n");
	
	/* Volume vector */
    printf("vol:");
    for(i=0;i<model->Ncells;i++){
        printf("%f ",model->vol[i]);
    }
    printf("\n");

    
    for(i=0;i<model->Mspecies;i++){
        for (j=model->jcN[i]; j<model->jcN[i+1]; j++){
            printf("N: %i\t%lu - %lu\t%lu\t%i\n",j,(long unsigned)model->jcN[i],(long unsigned)model->jcN[i+1],(long unsigned)model->irN[j],model->prN[j]);
        }
    }
    printf("end N\n");

    for(i=0;i<model->Mreactions;i++){
        for (j=model->jcG[i]; j<model->jcG[i+1]; j++){
            printf("G %i\t%lu\n",j,(long unsigned)model->irG[j]);
        }
    }
    printf("end G\n");
	
	/* Initial condition */
    printf("u0:");
    for(i=0;i<Ndofs;i++){
        printf("%i ",model->u0[i]);
    }
    printf("\nend u0\n");
	
	/* Time span 
	int tlen;
	double *tspan;*/
    printf("tspan(%i):",model->tlen);
    for(i=0;i<model->tlen;i++){
        printf("%f ",model->tspan[i]);
    }
    printf("\nend tspan\n");
	
	/* Subdomain vector 
	int *sd;*/
    printf("sd:");
    for(i=0;i<model->Ncells;i++){
        printf("%i ",model->sd[i]);
    }
    printf("\nend sd\n");
	
	/* Data vector 
	int dsize;
	double *data;*/
    printf("data:");
    for(i=0;i<model->Ncells*model->dsize;i++){
        printf("%f ",model->data[i]);
    }
    printf("\nend data\n");
    //------------
}
