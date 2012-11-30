 /* Propensity definition for the 'brianmike' system. */

#include <stdlib.h>
#include <stdio.h>
/* Type definition for the propensity functions. */
#include "propensities.h"

/* Number of reactions */
#define NR	8

const double N_A = 6.022e23;
const double SA = 201.056; /* CHECK!! surface area scaling ?? */

/* Rate constants */ 
const double kRs  = 4.0;  
const double kRd0 = 4e-4;
const double kRL  = 2e-3;
const double kRLm = 1e-2;
const double kRd1 = 4e-4;
const double kGa  = 1e-6; 
const double kGd  = 0.1; 
const double kG1  = 1.0;  

// Species
// R RL G Ga Gbg Gd
#define s_R   0
#define s_RL  1
#define s_G   2
#define s_Ga  3
#define s_Gbg 4
#define s_Gd  5

/* Reaction propensities. */
double rFun1(const int *x, double t, double vol, const double *data, int sd)
/* EmptySet -> R */
{
  return kRs/SA*vol;
}

double rFun2(const int *x, double t, double vol, const double *data, int sd)
/* R -> EmptySet */
{
  return kRd0*x[s_R];
}

double rFun3(const int *x,  double t,double vol, const double *data, int sd)
/* L + R -> RL + L. Fix the gradient expression using the data vector. */
{
    //#define MOLAR  (6.02e-01*vol)
	//return kRL/MOLAR*x[s_R]*data[0]/vol;
	return kRL*x[s_R]*data[0]/vol/10.0;
}

double rFun4(const int *x, double t, double vol, const double *data, int sd)
/*RL -> R */
{
  return kRLm*x[s_RL];
}

double rFun5(const int *x, double t, double vol, const double *data, int sd)
/*RL -> R */
{
  return kRd1*x[s_RL];
}

double rFun6(const int *x, double t, double vol, const double *data, int sd)
/*RL + G -> Ga + Gbg + RL*/
{
  return SA*kGa*x[s_RL]*x[s_G]/vol; 
}

double rFun7(const int *x, double t, double vol, const double *data, int sd)
/*Ga -> Gd*/
{
  return kGd*x[s_Ga];
}

double rFun8(const int *x, double t, double vol, const double *data, int sd)
/*Gd + Gbg -> G*/
{
  return SA*kG1*x[s_Gd]*x[s_Gbg]/vol;
}

PropensityFun *ALLOC_propensities(void)
{
  PropensityFun *ptr = malloc(sizeof(PropensityFun)*NR);
  
  ptr[0]=rFun1;
  ptr[1]=rFun2;
  ptr[2]=rFun3;
  ptr[3]=rFun4;
  ptr[4]=rFun5;
  ptr[5]=rFun6;
  ptr[6]=rFun7;
  ptr[7]=rFun8;
 
  
  return ptr;
}

void FREE_propensities(PropensityFun* ptr)
{
  free(ptr);
}


