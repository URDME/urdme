#ifndef DFSP__H
#define DFSP__H
#include "propensities.h"
#include "report.h"


void dfsp_core(const size_t *irD, const size_t *jcD, const double *prD,
            const int *u0, const size_t *irN, const size_t *jcN, const int *prN,
            const size_t *irG, const size_t *jcG, const double *tspan, const size_t tlen, int *U,
            const double *vol, const double *data, const int *sd, const size_t Ncells,
            const size_t Mspecies, const size_t Mreactions, const size_t dsize,
            int report_level, const double tau_d);

int dfsp_diffusion(int *xx, const size_t *irD, const size_t *jcD, const double *prD,
            const size_t Ncells,int Mspecies);

int dfsp_reactions(int *xx,
            const size_t *irN, const size_t *jcN, const int *prN,
            const size_t *irG, const size_t *jcG,
            PropensityFun *rfun,
            double tau_d, double tt,
            const double *vol, const double *data, const int *sd,
            const size_t Ncells,int Mspecies,int Mreactions,const size_t dsize);

#endif
