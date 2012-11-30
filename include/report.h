/* report.h -- Typedef for report functions to be used with nsm_core. */
/* A. Hellander 2008-12-08 */

#ifndef REPORT__H
#define REPORT__H

/* Report function datatype */ 
typedef void (*ReportFun)(double, const double, const double, long int, long int, int, int);

/* Provided report function */
void reportFun1(double, const double, const double, long int,long int, int, int);

#endif /* REPORT__H */
