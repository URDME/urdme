/* spin.c */

/* P. Bauer   2012-04-02 */
/* J. Cullhed 2008-08-04. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mcheck.h>
#include "../src/nsm/binheap.h"
#include "../include/propensities.h"
#include "../src/nsm/nsm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define N_TESTS 10

/* External variable used to select correct propensity function */
int testNumber;

/*---------------------------------------------------------------------*/
int is_heap(const double *data, const int *i1, const int *i2,
    const size_t len, const double *reaction_table)
/* Check if data, i1 and i2 constructs a valid binary heap. */
{
  size_t i;
  /* Test structure on data. */
  for (i=1; i<(len-1)/2; i++)
    if (data[i-1]>data[2*i-1] || data[i-1]>data[2*i]) return 0;

  /* Test structure on i1. */
  for (i=0; i<len; i++)
    if (reaction_table[i1[i]]!=data[i]) return 0;

  /* Test structure on i2. */
  for (i=0; i<len; i++)
    if (reaction_table[i]!=data[i2[i]]) return 0;

  return 1;
}
/*---------------------------------------------------------------------*/
int spin1(void)
/* Test of reaction-diffusion in four cells with one reaction: rFun1. */
{
  int ok=1;

  /* Length constants. */
  const size_t Ncells=4;
  const size_t Mspecies=3;
  const size_t Mreactions=1;
  const size_t Ndofs=Ncells*Mspecies;
  const size_t Dim=2;
  size_t i; 
  //size_t j;

  int xsum=0, ysum=0, zsum=0;

  /* Diffusion matrix D in sparse format. */
  const size_t irD[]={0,3,6,1,4,7,2,5,8,
                   0,3,9,1,4,10,2,5,11,
                   0,6,9,1,7,10,2,8,11,
                   3,6,9,4,7,10,5,8,11};

  const size_t jcD[]={0,3,6,9,12,15,18,21,24,27,30,33,36};

  const double prD[]={-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,
                      1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,
                      1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,
                      1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0};

  /* Initial vector u0. */
  const int u0[]={3,3,0,3,3,0,3,3,0,2,1,0};

  /* Stochiometric matrix N in sparse format. */
  const size_t irN[]={0,1,2};
  const size_t jcN[]={0,3};
  const int prN[]={-1,-1,1};

  /* Dependency graph G in sparse format. */
  const size_t irG[]={0,0,0};
  const size_t jcG[]={0,1,2,2,3};

  /* Output times tspan. */
  const double tspan[]={0,1,2,30,40,53,61,70,100,109,231};
  const size_t tlen=sizeof(tspan)/sizeof(tspan[0]);

  /* Output U. */
  int *U=(int *)malloc(Ndofs*tlen*sizeof(int));

  /* Other input. vol, pos and sd. */
  const double vol[]={1.0,1.0,1.0,1.0};
  const int sd[]={1,1,1,1};
  
  /* Propensity. */
  testNumber=1;
  mtrace();
  nsm_core(irD,jcD,prD,u0,irN,jcN,prN,irG,jcG,tspan,tlen,U,vol,NULL,sd,Ncells,
          Mspecies,Mreactions,Dim,0);
  muntrace();

  /* Check result. */
  for (i=0; i<Ndofs; i+=Mspecies) xsum+=U[(tlen-1)*Ndofs+i];
  for (i=1; i<Ndofs; i+=Mspecies) ysum+=U[(tlen-1)*Ndofs+i];
  for (i=2; i<Ndofs; i+=Mspecies) zsum+=U[(tlen-1)*Ndofs+i];

  if (xsum!=1 || ysum!=0 || zsum!= 10) ok=0;

  /* Always check for negative states. */
  for (i=0; i<Ndofs*tlen; i++) ok=ok && U[i]>=0;

  /*
  for (i=0; i<tlen; i++) {
    for (j=0; j<Ndofs; j+=3)
      printf("%i %i %i   ", U[i*Ndofs+j], U[i*Ndofs+j+1], U[i*Ndofs+j+2]);
    printf("\n");
  }
  */

  free(U);

  return ok;
}
/*---------------------------------------------------------------------*/
int spin2(void)
/* Test of reaction-diffusion in four cells with two reactions: rFun2 and 
   rFun3. */
{
  int ok=1;

  /* Length constants. */
  const size_t Ncells=4;
  const size_t Mspecies=3;
  const size_t Mreactions=2;
  const size_t Ndofs=Ncells*Mspecies;
  const size_t Dim=2;
  size_t i;
  //size_t j;

  int xsum=0, ysum=0, zsum=0;

  /* Diffusion matrix D in sparse format. */
  const size_t irD[]={0,3,6,1,4,7,2,5,8,
                   0,3,9,1,4,10,2,5,11,
                   0,6,9,1,7,10,2,8,11,
                   3,6,9,4,7,10,5,8,11};

  const size_t jcD[]={0,3,6,9,12,15,18,21,24,27,30,33,36};

  const double prD[]={-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,
                      1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,
                      1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,
                      1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0};

  /* Initial vector u0. */
  const int u0[]={10,5,0,30,20,0,5,15,0,5,10,0};

  /* Stochiometric matrix N in sparse format. */
  const size_t irN[]={0,2,1,2};
  const size_t jcN[]={0,2,4};
  const int prN[]={-2,1,-2,1};

  /* Dependency graph G in sparse format. */
  const size_t irG[]={0,1,0,1};
  const size_t jcG[]={0,1,2,2,3,4};

  /* Output times tspan. */
  const double tspan[]={0,1,2,30,40,53,61,70,100,109,231};
  const unsigned tlen=sizeof(tspan)/sizeof(tspan[0]);

  /* Output U. */
  int *U=(int *)malloc(Ndofs*tlen*sizeof(int));

  /* Other input. vol, pos and sd. */
  const double vol[]={1.0,1.0,1.0,1.0};
  const int sd[]={1,1,1,1};

  /* Propensity. */
  testNumber=2;

  nsm_core(irD,jcD,prD,u0,irN,jcN,prN,irG,jcG,tspan,tlen,U,vol,NULL,sd,Ncells,
          Mspecies,Mreactions,Dim,0);

  /* Check result. */
  for (i=0; i<Ndofs; i+=Mspecies) xsum+=U[(tlen-1)*Ndofs+i];
  for (i=1; i<Ndofs; i+=Mspecies) ysum+=U[(tlen-1)*Ndofs+i];
  for (i=2; i<Ndofs; i+=Mspecies) zsum+=U[(tlen-1)*Ndofs+i];

  if (xsum!=0 || ysum!=0 || zsum!= 50) ok=0;

  /* Always check for negative states. */
  for (i=0; i<Ndofs*tlen; i++) ok=ok && U[i]>=0;

  /*
  for (i=0; i<tlen; i++) {
    for (j=0; j<Ndofs; j+=3)
      printf("%i %i %i   ", U[i*Ndofs+j], U[i*Ndofs+j+1], U[i*Ndofs+j+2]);
    printf("\n");
  }
  */

  free(U);

  return ok;
}
/*---------------------------------------------------------------------*/
int spin3(void)
/* Test of reaction-diffusion in four cells with four reactions: rFun2, rFun3
   and rFun4. */
{
  int ok=1;

  /* Length constants. */
  const size_t Ncells=4;
  const size_t Mspecies=3;
  const size_t Mreactions=3;
  const size_t Ndofs=Ncells*Mspecies;
  const size_t Dim=2;
  size_t i;
  //size_t j;

  int xsum=0, ysum=0, zsum=0;

  /* Diffusion matrix D in sparse format. */
  const size_t irD[]={0,3,6,1,4,7,2,5,8,
                   0,3,9,1,4,10,2,5,11,
                   0,6,9,1,7,10,2,8,11,
                   3,6,9,4,7,10,5,8,11};

  const size_t jcD[]={0,3,6,9,12,15,18,21,24,27,30,33,36};

  const double prD[]={-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,
                      1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,
                      1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,
                      1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0};

  /* Initial vector u0. */
  const int u0[]={10,25,0,10,5,0,10,5,0,20,15,0};

  /* Stochiometric matrix N in sparse format. */
  const size_t irN[]={0,2,1,2,2};
  const size_t jcN[]={0,2,4,5};
  const int prN[]={-2,1,-2,1,-1};

  /* Dependency graph G in sparse format. */
  const size_t irG[]={0,1,2,0,2,1,2,2};
  const size_t jcG[]={0,1,2,3,5,7,8};

  /* Output times tspan. */
  const double tspan[]={0,1,2,30,40,53,61,70,100,109,231};
  const unsigned tlen=sizeof(tspan)/sizeof(tspan[0]);

  /* Output U. */
  int *U=(int *)malloc(Ndofs*tlen*sizeof(int));

  /* Other input. vol, pos and sd. */
  const double vol[]={1.0,1.0,1.0,1.0};
  const int sd[]={1,1,1,1};

  /* Propensity. */
  testNumber=3;
  
  nsm_core(irD,jcD,prD,u0,irN,jcN,prN,irG,jcG,tspan,tlen,U,vol,NULL,sd,Ncells,
          Mspecies,Mreactions,Dim,0);  

  /* Check result. */
  for (i=0; i<Ndofs; i+=Mspecies) xsum+=U[(tlen-1)*Ndofs+i];
  for (i=1; i<Ndofs; i+=Mspecies) ysum+=U[(tlen-1)*Ndofs+i];
  for (i=2; i<Ndofs; i+=Mspecies) zsum+=U[(tlen-1)*Ndofs+i];

  if (xsum!=0 || ysum!=0 || zsum!= 0) ok=0;

  /* Always check for negative states. */
  for (i=0; i<Ndofs*tlen; i++) ok=ok && U[i]>=0;

  /*
  for (i=0; i<tlen; i++) {
    for (j=0; j<Ndofs; j+=3)
      printf("%i %i %i   ", U[i*Ndofs+j], U[i*Ndofs+j+1], U[i*Ndofs+j+2]);
    printf("\n");
  }
  */

  free(U);

  return ok;
}
/*---------------------------------------------------------------------*/
int spin4(void)
/* Test of reaction-diffusion in four cells with one reaction: rFun5.
   With two subdomains. No diffusion from cell zero. */
{
  int ok=1;

  /* Length constants. */
  const size_t Ncells=4;
  const size_t Mspecies=3;
  const size_t Mreactions=2;
  const size_t Ndofs=Ncells*Mspecies;
  const size_t Dim=2;
  size_t i;
  //size_t j;

  //int xsum=0, ysum=0, zsum=0;

  /* Diffusion matrix D in sparse format. */
  const size_t irD[]={0,1,2,0,3,9,1,4,10,2,5,11,
                    0,6,9,1,7,10,2,8,11,
                    3,6,9,4,7,10,5,8,11};

  const size_t jcD[]={0,1,2,3,6,9,12,15,18,21,24,27,30};

  const double prD[]={-2.0,-2.0,-2.0,
                      1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,
                      1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,
                      1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0};

  /* Initial vector u0. */
  const int u0[]={1,1,1,3,3,3,1,2,3,5,4,3};

  /* Stochiometric matrix N in sparse format. */
  const size_t irN[]={0,1,2,0,1,2};
  const size_t jcN[]={0,3,6};
  const int prN[]={-1,-1,1,1,1,-1};

  /* Dependency graph G in sparse format. */
  const size_t irG[]={0,0,1,0,1,0,1};
  const size_t jcG[]={0,1,2,3,5,7};

  /* Output times tspan. */
  const double tspan[]={0,1,2,30,40,53,61,70,100,109,231};
  const size_t tlen=sizeof(tspan)/sizeof(tspan[0]);

  /* Output U. */
  int *U=(int *)malloc(Ndofs*tlen*sizeof(int));

  /* Other input. vol, pos and sd. */
  const double vol[]={1.0,1.0,1.0,1.0};
  const int sd[]={1,0,0,0};

  /* Propensity. */
  testNumber=4;
  
  nsm_core(irD,jcD,prD,u0,irN,jcN,prN,irG,jcG,tspan,tlen,U,vol,NULL,sd,Ncells,
          Mspecies,Mreactions,Dim,0);  

  /* Check result. */
  if(U[(tlen-1)*Ndofs+2]!=20) ok=0;

  /* Always check for negative states. */
  for (i=0; i<Ndofs*tlen; i++) ok=ok && U[i]>=0;

  /*
  for (i=0; i<tlen; i++) {
    for (j=0; j<Ndofs; j+=3)
      printf("%i %i %i   ", U[i*Ndofs+j], U[i*Ndofs+j+1], U[i*Ndofs+j+2]);
    printf("\n");
  }
  */

  free(U);

  return ok;
}
/*---------------------------------------------------------------------*/
int spin5(void)
/* Test check if negative number of molecules. */
{
  int ok=1;

  /* Length constants. */
  const size_t Ncells=4;
  const size_t Mspecies=3;
  const size_t Mreactions=1;
  const size_t Ndofs=Ncells*Mspecies;
  const size_t Dim=2;
  size_t i, j;

  //int xsum=0, ysum=0, zsum=0;

  /* Diffusion matrix D in sparse format. */
  const size_t irD[]={0,3,6,1,4,7,2,5,8,
                   0,3,9,1,4,10,2,5,11,
                   0,6,9,1,7,10,2,8,11,
                   3,6,9,4,7,10,5,8,11};

  const size_t jcD[]={0,3,6,9,12,15,18,21,24,27,30,33,36};

  const double prD[]={-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,
                      1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,
                      1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,
                      1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0};

  /* Initial vector u0. */
  const int u0[]={10,10,0,0,0,0,0,0,0,0,0,0};

  /* Stochiometric matrix N in sparse format. */
  const size_t irN[]={0,1,2};
  const size_t jcN[]={0,3};
  const int prN[]={-1,-1,1};

  /* Dependency graph G in sparse format. */
  const size_t irG[]={0,0,0};
  const size_t jcG[]={0,1,2,2,3};

  /* Output times tspan. */
  const double tspan[]={0,1,2,30,40,53,61,70,100,109,231};
  const size_t tlen=sizeof(tspan)/sizeof(tspan[0]);

  /* Output U. */
  int *U=(int *)malloc(Ndofs*tlen*sizeof(int));

  /* Other input. vol, pos and sd. */
  const double vol[]={1.0,1.0,1.0,1.0};
  const int sd[]={1,1,1,1};

  /* Inline propensities. */
/* const double K[]={0.0,0.0,1.0};
  const int I[]={0,0,0};
  const unsigned jcS[]={0,1};
  const int prS[]={0}; */
  testNumber=5;
  
  nsm_core(irD,jcD,prD,u0,irN,jcN,prN,irG,jcG,tspan,tlen,U,vol,NULL,sd,Ncells,
          Mspecies,Mreactions,Dim,0);  

  int neg=0;
  for (i=0; i<tlen; i++) {
    for (j=0; j<Ndofs; j++) {
      if (neg)
        ok=ok && U[(tlen-1)*Ndofs+i]==0;
      else if (U[i*Ndofs+j]<0) {
        neg=1;
        break;
      }
    }
  }
  free(U);
  return ok;
}
/*---------------------------------------------------------------------*/
int spin6(void)
/* Binary heap spin. Test initialize_heap. */
{
  double *times, *reaction_table;
  int *i1, *i2;
  size_t i, len=1;
  int ok=1;

  while (len<1000) {
    times=(double *)malloc(len*sizeof(double));
    reaction_table=(double *)malloc(len*sizeof(double));
    i1=(int *)malloc(len*sizeof(int));
    i2=(int *)malloc(len*sizeof(int));

    for (i=0; i<len; i++) {
      times[i]=drand48();
      reaction_table[i]=times[i];
      i1[i]=i;
      i2[i]=i;
    }

    initialize_heap(times, i1, i2, len);

    ok=ok && is_heap(times, i1, i2, len, reaction_table);

    free(i2);
    free(i1);
    free(times);
    len++;
  }
  return ok;
}
/*---------------------------------------------------------------------*/
int spin7(void)
/* Binary heap spin. Test update. */
{
  double *times, *reaction_table;
  int *i1, *i2;
  size_t i, len=1;
  int ok=1;

  while (len<1000) {
    times=(double *)malloc(len*sizeof(double));
    reaction_table=(double *)malloc(len*sizeof(double));
    i1=(int *)malloc(len*sizeof(int));
    i2=(int *)malloc(len*sizeof(int));

    for (i=0; i<len; i++) {
      times[i]=drand48();
      reaction_table[i]=times[i];
      i1[i]=i;
      i2[i]=i;
    }

    initialize_heap(times, i1, i2, len);

    /* Update heap with random values. */
    for (i=0; i<len; i++) {
      const double t=drand48();
      times[i2[i]]=t;
      reaction_table[i]=t;
      update(i2[i], times, i1, i2, len);
    }

    ok=ok && is_heap(times, i1, i2, len, reaction_table);

    free(i2);
    free(i1);
    free(times);
    len++;
  }
  return ok;
}
/*---------------------------------------------------------------------*/
int spin8(void)
/* Binary heap spin. Test percolate down, using update.*/
{
  double *times, *reaction_table;
  int *i1, *i2;
  size_t i, len=100;
  int ok=1;

  times=(double *)malloc(len*sizeof(double));
  reaction_table=(double *)malloc(len*sizeof(double));
  i1=(int *)malloc(len*sizeof(int));
  i2=(int *)malloc(len*sizeof(int));

  for (i=0; i<len; i++) {
    times[i]=drand48();
    reaction_table[i]=times[i];
    i1[i]=i;
    i2[i]=i;
  }

  initialize_heap(times, i1, i2, len);

  /* Insert one big number that should be in the lowest level. */
  times[i2[0]]=10.0;
  update(i2[0], times, i1, i2, len);
  reaction_table[0]=10.0;

  ok=ok && is_heap(times, i1, i2, len, reaction_table);

  free(i2);
  free(i1);
  free(times);
  len++;
  return ok;
}
/*---------------------------------------------------------------------*/
int spin9(void)
/* Binary heap spin. Test percolate up, using update.*/
{
  double *times, *reaction_table;
  int *i1, *i2;
  size_t i, len=100;
  int ok=1;

  times=(double *)malloc(len*sizeof(double));
  reaction_table=(double *)malloc(len*sizeof(double));
  i1=(int *)malloc(len*sizeof(int));
  i2=(int *)malloc(len*sizeof(int));

  for (i=0; i<len; i++) {
    times[i]=drand48();
    reaction_table[i]=times[i];
    i1[i]=i;
    i2[i]=i;
  }

  initialize_heap(times, i1, i2, len);

  /* Insert one small number that should be in the top node. */
  times[i2[len-1]]=-10.0;
  reaction_table[len-1]=-10.0;
  update(i2[len-1], times, i1, i2, len);

  ok=ok && is_heap(times, i1, i2, len, reaction_table);

  /* Check that -10.0 is in the top node. */
  ok=ok && (times[0]==-10.0);

  free(i2);
  free(i1);
  free(times);
  len++;
  return ok;
}
/*---------------------------------------------------------------------*/
int spin10(void)
/* Binary heap spin. Messy update test. */
{
  double *times, *reaction_table;
  int *i1, *i2;
  size_t i, len=1000;
  int ok=1;

  times=(double *)malloc(len*sizeof(double));
  reaction_table=(double *)malloc(len*sizeof(double));
  i1=(int *)malloc(len*sizeof(int));
  i2=(int *)malloc(len*sizeof(int));

  for (i=0; i<len; i++) {
    times[i]=drand48();
    reaction_table[i]=times[i];
    i1[i]=i;
    i2[i]=i;
  }

  initialize_heap(times, i1, i2, len);

   /* Update heap with random values. */
  for (i=0; i<len; i++) {
    const double t=drand48();
    times[i2[i]]=t;
    reaction_table[i]=t;
    update(i2[i], times, i1, i2, len);
  }

  /* Update heap with random values at random places. */
  for (i=0; i<len; i++) {
    const double t=drand48();
    const size_t pos=rand()%len;
    times[i2[pos]]=t;
    reaction_table[pos]=t;
    update(i2[pos], times, i1, i2, len);
  }

  ok=ok && is_heap(times, i1, i2, len, reaction_table);

  free(i2);
  free(i1);
  free(times);
  len++;
  return ok;
}
/*---------------------------------------------------------------------*/
int main(void) {
  srand48(1234);
  size_t i;
  int (*f[N_TESTS])()={spin1, spin2, spin3, spin4, spin5, spin6, spin7,spin8, spin9, spin10};
  int passed=0;

  /* Run all specified tests. */
  for (i=0; i<N_TESTS; i++) {
    if (f[i]()) {
      printf("Spin #%zu passed.\n", i+1); 
      passed++;
    } else {
      printf("Spin #%zu failed.\n", i+1);
    }
  }
  printf("Passed: %d/%d test(s).\n",passed, N_TESTS);
}
/*---------------------------------------------------------------------*/

