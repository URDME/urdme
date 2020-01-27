/* NMDAfun.h - Auxiliary function declarations for NMDA synaptic propensities. */

/* S. Engblom 2019-11-29 */

#ifndef NMDAFUN_H
#define NMDAFUN_H

const static double valence = -2.0;
const static double memb_fraction = 0.8;

/* forward declarations */
double rmb(const double v);
double rmu(const double v);

#endif /* NMDAFUN_H */
