#ifndef _LOGLIKE_PSI_H_
#define _LOGLIKE_PSI_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



void 	loglike_psi_initialize();
void 	loglike_psi_finalize();
float 	loglike_psi(float *x, int n, void *output, int *info);

float log_priorHB(float *x, float *psi, int n);


#endif
