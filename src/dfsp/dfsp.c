/* Stand-alone dfsp solver for use with URDME. */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "mat.h"
#include "mex.h"
#include "propensities.h"
#include "dfsp.h"
#include "matmodel.h"


void *dfsp(void *data);

int main(int argc, char *argv[])

{
	
	char *infile,*outfile;
	int i, nt=1, report_level=1;
	
	/* TODO. Parsing of input with error checking. 
	   Input syntax on form: -numtraj=4 etc. ?*/
	
	if (argc < 3){
		printf("To few arguments to nsm.");
	    exit(-1);	
	}
	
	/* Input file. */
	infile  = argv[1];
	
	/* Output file (or directory). */
	outfile = argv[2]; 
	
	/* Number of trajectories (realizations) to generate. This 
	   feature is disabled in the main Matlab interface 
	   in this release. */
	
	if (argc > 3)
		nt = atoi(argv[3]);
	
	/* Report level (passed to the core solver). 
	   Defaults to 1 (partial report) for single trajectory runs,
	   and 0 (no report) for multi-trajectory runs. */
	if (argc > 4)
		report_level = atoi(argv[4]);
	else if (nt > 1)
		report_level = 0;
	
	
	/* Read model specification */
	urdme_model *model;
	model = read_model(infile);
	model->infile = infile;
	
	if (model == NULL){
		printf("Fatal error. Couldn't load model file or currupt model file.");
		return(-1);
	}
	
	/* Initialize RNG(s).  */
    srand48( time(NULL) * getpid() );
	

    /* Initialize extra args */
	model->num_extra_args = 5;
	model->extra_args=malloc(model->num_extra_args*sizeof(void *));
	
	/* Set report level */
	model->extra_args[0] = malloc(sizeof(int));
	*(int *)(model->extra_args[0]) = report_level;

    /*reopen MAT file*/
    MATFile *input_file;
    mxArray *DT,*sopts;

    input_file = matOpen(infile,"r");
    if (input_file == NULL){
        printf("Failed to open mat-file.\n");
        return(-1);   
    }
    /* Set tau_d */
    sopts = matGetVariable(input_file, "sopts");
    if(sopts==NULL){
        printf("Step length (tau_d) is missing in the model file\n");
        return(-1);   
    }
    model->extra_args[1] = malloc(sizeof(double));
    *(double *)(model->extra_args[1]) = *(double *)mxGetPr(sopts);
    //fprintf(stderr,"here tau_d=%e\n",*(double *)(model->extra_args[1]));

    /* Set DT, DFSP transition matrix */
    DT = matGetVariable(input_file, "DT");
    if (DT == NULL){
        printf("Step length (tau_d) is missing in the model file\n");
        return(-1);   
    }
    model->extra_args[2]  = (size_t *)mxGetIr(DT);
    model->extra_args[3]  = (size_t *)mxGetJc(DT);
    model->extra_args[4]  = (double *)mxGetPr(DT);
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
	
	//destroy_model(model);
	
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


 
