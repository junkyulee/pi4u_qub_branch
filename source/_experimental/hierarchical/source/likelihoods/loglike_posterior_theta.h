#ifndef _LOGLIKE_POSTERIOR_THETA_FAST_H_
#define _LOGLIKE_POSTERIOR_THETA_FAST_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void 	loglike_posterior_theta_initialize();
void 	loglike_posterior_theta_finalize();
float 	loglike_posterior_theta(float *x, int n, void *output, int *info );

float 	loglike_(float *x, int n, void *output, int *info );

#endif
