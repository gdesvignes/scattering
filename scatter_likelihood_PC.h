//#ifndef SCATTERING_PC_H_INCLUDED
//#define SCATTERING_PC_H_INCLUDED

#include "scattering.h"

extern MNStruct *sp;

double ScatterLike_PC (double theta[], int nDims, double phi[], int nDerived, void *context);
void prior (double cube[], double theta[], int nDims, void * context);
void setup_loglikelihood();

//#else
//extern MNStruct *sp;

//#endif

