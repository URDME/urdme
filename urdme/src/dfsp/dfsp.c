/* Stand-alone dfsp solver for use with URDME. */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "propensities.h"
#include "dfsp.h"
#include "matmodel.h"

#ifndef URDME_LIBMAT
#include "mat.h"
#include "mex.h"
#else
#include "read_matfile.h"
#endif


void *dfsp(void *data);

int main(int argc, char *argv[])

{
	
	char *infile,*outfile;
	int i, nt=1, report_level=1;
	
	if (argc < 3){
		printf("To few arguments to dfsp.");
	    exit(-1);	
	}
	
	/* Input file. */
	infile  = argv[1];
	
	/* Output file (or directory). */
	outfile = argv[2]; 

	/* Read model specification */
	urdme_model *model;
	model = read_model(infile);
	model->infile = infile;

	if (model == NULL){
		printf("Fatal error. Couldn't load model file or currupt model file.");
		return(-1);
	}
	
    /* Initialize extra args */
	model->num_extra_args = 5;
	model->extra_args=malloc(model->num_extra_args*sizeof(void *));

    /*reopen MAT file*/
    MATFile *input_file;
    mxArray *DT,*sopts;
    
    input_file = matOpen(infile,"r");
    if (input_file == NULL){
        printf("Failed to open mat-file.\n");
        return(-1);
    }
    
    /* Look for seed */
    mxArray *mxseed;
	mxseed = matGetVariable(input_file, "seed");
    if(mxseed!=NULL && (!mxIsEmpty(mxseed))) {
        srand48((long int)mxGetScalar(mxseed));
    } else {
        srand48((long int)time(NULL)+(long int)(1e9*clock()));
    }

	/*
     If seed is provided as a parameter, it takes precedence.
     We need to be able to pass the seed as a paramters as well
     as in the input file in the cases where the solver is run
     in a distributed environment.
     */
	if (argc > 3) {
		srand48((long int)atoi(argv[3]));
	}

    /* Look for an optional parameter matrix. */
	const double *matfile_parameters;
	int mpar = 0;
	int npar = 0;
	int param_case=1;
    mxArray *mxparameters;
	mxparameters = matGetVariable(input_file, "parameters");
    if (mxparameters != NULL) {
		matfile_parameters = (double *)mxGetPr(mxparameters);
		mpar = mxGetM(mxparameters);
		npar = mxGetN(mxparameters);
	}

    /* Look if a parameter case if supplied as a parameter. */
	if (argc > 4) {
	    param_case = (int)atoi(argv[4]);
	}
	
	if (param_case > npar && mxparameters!=NULL) {
		perror("dfspcore: Fatal error, parameter case is larger than n-dimension in parameter matrix.\n");
		exit(-2);
	}
	
	/* Create global parameter variable for this parameter case. */
    parameters = (double *)malloc(mpar*sizeof(double));
    memcpy(parameters,&matfile_parameters[mpar*(param_case-1)],mpar*sizeof(double));
    
    mxArray *mxreport;
 
	mxreport = matGetVariable(input_file, "report");
	if (mxreport != NULL)
		report_level = (int) mxGetScalar(mxreport);
	else
        report_level=1;
    
    /* Set report level */
	model->extra_args[0] = (int *)malloc(sizeof(int));
	*(int *)(model->extra_args[0]) = report_level;
    
    /* Set tau_d */
    sopts = matGetVariable(input_file, "sopts");
    if(sopts==NULL||mxIsEmpty(sopts)){
        printf("Step length (tau_d) is missing in the model file\n");
        return(-1);   
    }
    model->extra_args[1] = malloc(sizeof(double));
    *(double *)(model->extra_args[1]) = *(double *)mxGetPr(sopts);

    /* Set DT, DFSP transition matrix */
    DT = matGetVariable(input_file, "DT");
    if (DT == NULL){
        printf("Step length (tau_d) is missing in the model file\n");
        return(-1);   
    }
	
    /* Since the extra arg is freed in destroy model, we need to create 
       copies of these arrays to make sure that all memory in the model struct
       is allocated with malloc, not mxMalloc. */
    int Ndofs = model->Mspecies*model->Ncells;
    int nnzDT = mxGetNzmax(DT);
    size_t *irDT,*jcDT;
    double *prDT; 
    irDT = (size_t *)malloc(nnzDT*sizeof(size_t));	 
    jcDT = (size_t *)malloc((Ndofs+1)*sizeof(size_t));	
    prDT = (double *)malloc(nnzDT*sizeof(double));	
    memcpy(irDT,(size_t *)mxGetIr(DT),nnzDT*sizeof(size_t));
    memcpy(jcDT,(size_t *)mxGetJc(DT),(Ndofs+1)*sizeof(size_t));
    memcpy(prDT,(double *)mxGetPr(DT),nnzDT*sizeof(double));
    model->extra_args[2]  = irDT;
    model->extra_args[3]  = jcDT;
    model->extra_args[4]  = prDT;


    /* close MAT file*/
    matClose(input_file);

	/* Allocate memory to hold nt solutions. */
	init_sol(model,nt);
	
	/* Call nsm-solver: get a trajectories (will be threaded in a near future). */
	for (i=0;i<nt;i++){
		dfsp(model);
	}
	
	/* Print result to file(s)*/
	int didprint;
	didprint=dump_results(model, outfile,"single");
	if (!didprint){
		printf("Error: Failed to print results to file.\n");
		exit(-1);
	}
	
    free(parameters);
	destroy_model(model);
	return(0);
	
}


void *dfsp(void *data){
	
	/* Unpack input */
	urdme_model* model;
	model = (urdme_model *)data;
	int Ndofs, *U;

    /* sopts(1) = tau_d */
    double tau_d = *(double *)model->extra_args[1];
    if(tau_d < 0 || tau_d > model->tspan[1]){
        printf("Warning, tau_d is not properly defined, using %e\n",model->tspan[1]);
        tau_d = model->tspan[1];
    }

	/* Uses a report function with optional report level. This is
	 passed as the first extra argument. */ 
	int report_level = *(int *)model->extra_args[0];
	
	/* Output array (to hold a single trajectory) */
	Ndofs = model->Ncells*model->Mspecies;
	U = malloc(model->tlen*Ndofs*sizeof(int));
	
	/* Call the core simulation routine.  */
	dfsp_core((size_t *)model->extra_args[2], (size_t *)model->extra_args[3], (double *)model->extra_args[4],
             model->u0,
			 model->irN, model->jcN, model->prN, model->irG,
			 model->jcG, model->tspan, model->tlen, 
			 U, model->vol, model->data, model->sd, model->Ncells, 
			 model->Mspecies, model->Mreactions, model->dsize,
             report_level, tau_d);
	
	/* Attach trajectory to model. */
	if (model->nsol >= model->nsolmax){
		printf("Warning, you are trying to add more trajectories to model struct\n" );
	    printf("but it already contains its maximum number. This is not my fault, you screwed up.");
		printf("I'm starting to overwrite old ones (starting with the first one).\n");
		model->nsol = 0;
	}
	
	if ((model->U[model->nsol]) != NULL)
		free(model->U[model->nsol]);
	
	(model->U)[(model->nsol)++] = U;
	
	return U;
	
}


 
