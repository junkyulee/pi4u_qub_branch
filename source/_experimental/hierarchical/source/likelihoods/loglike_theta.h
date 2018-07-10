#ifndef _LOGLIKE_THETA_H_
#define _LOGLIKE_THETA_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void loglike_theta_initialize();
void loglike_theta_finalize();
float loglike_theta(float *x, int n, void *output, int *winfo);

float loglike_(float *x, int n, void *output, int *winfo);

#endif
