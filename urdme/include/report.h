/* report.h -- URDME report function. */

/* S. Engblom 2017-02-17 (Major revision, URDME 1.3, Comsol 5) */
/* A. Hellander 2008-12-08 */

#ifndef REPORT__H
#define REPORT__H

#include "mex.h"

#define PRINTF mexPrintf
#define PERROR mexErrMsgTxt

/* Report function datatype */ 
typedef void (*ReportFun)(double,double,double,long,long,int,int);

/* Provided report function */
void URDMEreportFun(double,double,double,long,long,int,int);

#endif /* REPORT__H */
