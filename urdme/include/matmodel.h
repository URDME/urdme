#ifndef MATMODEL__H
#define MATMODEL__H

/* 
   A. Hellander, 2010-04-20
   B. Drawert    2010-05-18
*/


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
 
 URDME model struct. 
 
 The fields of this struct is typically initialized by a call to "read_model". 
 The data structures are contained in the input file created by the 
 Matlab interface (.mat format)

 */

typedef struct{
	
	/* The "raw data" needed by the core nsm solver. 
	   These are essential fields, and they are initialized 
	   after a successful call to "read_model".  */
	
	/* Constants */
	int Mspecies;
	int Mreactions;
	int Ncells;

	/* Diffusion matrix (sparse CCS) */
	size_t * irD;
	size_t * jcD;
	double * prD;
	
	/* Volume vector */
	double *vol;
	
	/* Stoichiometry matrix (sparse CCS) */
	size_t* irN;
	size_t* jcN;
	int *prN;
	
	/* Dependency graph (sparse CCS) */
	size_t* irG;
	size_t* jcG;
	
	/* Initial condition */
	int *u0;
	
	/* Time span */
	int tlen;
	double *tspan;
	
	/* Subdomain vector */
	int *sd;
	
	/* Data vector */
	int dsize;
	double *data;
	
	/* 
     
       Output. It is up to the solver to attach the results after
	   a simulation to the urdme_struct. OBS! Solvers are neither
	   required nor guaranteed to do so, they may return a pointer 
	   to the result instead (or do both), if that's preferable for some reason.
     
     */
	
	int **U;
    
    
    
	/* Maximal number of solutions (set when U is allocated), defaults to one. */ 
	int nsolmax;
	/* Current number of stored solutions. */ 
	int nsol;
	
	/* 
	 
	   A solver might need additional data. Typically, this will be data that is not 
       directly related to the geometry and model (in constrast to the required fields above),
       but could be e.g. derived from those basic data types, or arguments that control the solver's behaviour. 
	   Examples could be the optional report level or a seed to the RNG that can be supplied to the NSM-solver, 
	   or the lookup-table and timestep used by the contributed plug-in DFSP solver.  
	 
	   OBS: this data can be stored in the input model file, 
       but "read_model" will not look for it and initialize it in the struct. 
	   Rather, your particular solver will need to know to look for it and extract it manually. 
	   
	   For example of use, see the core nsm solver (nsm.c) or the plug-in DFSP solver (dfsp.c). 
	 
	*/
	 
	int num_extra_args;
	void **extra_args;
	
	/* Name of input model file that the struct was generated from. */
	char *infile;
	
	
} urdme_model;

urdme_model *read_model(char *file);

void read_solution(urdme_model *model, char*file);
void init_sol(urdme_model *model, int nsolmax);

int destroy_model(urdme_model *model);
int dump_results(urdme_model* model, char *filename, char *type);

void debug_print_model(urdme_model* model);

#endif
