/* C-model file for the neuron example. */


#include <stdlib.h>
#include "propensities.h"

#define V   0
#define Vk  1
#define Vd  2


#define NR       8 	 

#define BULK     3 
#define MEMBRANE 2


/* Parameters. */

/* Number of V created per second in the whole domain. */
const double kb          = 10.0; 
const double kd          = 1.0; 
const double sigma_f     = 100.0; 
//const double sigma_kd    = 1.0; 
//const double sigma_dk    = 100.0; 
const double mu          = 0.1;


const double totvol = 6.17e-5; /* Total volume above z = 1.6 */
 
/* Reaction propensities. */


double rFun1(const int *x, double t, double vol, const double *data, int sd)
/* EmptySet -> V    
    
 * Creation of vescicle. 
    
 */

{
if ((int)data[0]==1) 
		return kb*vol/totvol; 
	return 0.0;
}

double rFun2(const int *x, double t, double vol, const double *data, int sd)
/* V -> Vk    
    
 * Binding of V to fiber (kinesin state).  
 * 
 */

{
   return sigma_f*x[V]; 
}

double rFun3(const int *x, double t, double vol, const double *data, int sd)
/* V -> Vd    
    
 * Binding of V to fiber (dynein state).  
 * 
 */

{
   return sigma_f*x[V]; 
}

double rFun4(const int *x, double t, double vol, const double *data, int sd)
/* Vk -> V    
    
 * Unbinding from fiber
  
 */

{
    return kd*x[Vk]; 
}

double rFun5(const int *x, double t, double vol, const double *data, int sd)
/* Vd -> V    
    
 * Unbinding from fiber
  
 */

{
    return kd*x[Vd]; 
}


/* Reversal of direction (active motor protein) */

double rFun6(const int *x, double t, double vol, const double *data, int sd)
/* Vk -> Vd*/
{
    double sigma_kd;
    
    if (t < 80.0)
        sigma_kd=1.0;
    else 
        sigma_kd=10.0;
    
	return sigma_kd*x[Vk];
}


double rFun7(const int *x, double t, double vol, const double *data, int sd)
/* Vk -> Vd*/
{
    double sigma_dk;
    if (t < 80.0)
        sigma_dk=10.0;
    else 
        sigma_dk=1.0;
    
	return sigma_dk*x[Vd];
}

double rFun8(const int *x, double t, double vol, const double *data, int sd)
/* V -> EmptySet*/
{
	return mu*x[V];
}



PropensityFun *ALLOC_propensities(void)
{
  PropensityFun *ptr = (PropensityFun *)malloc(sizeof(PropensityFun)*NR);
  
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



