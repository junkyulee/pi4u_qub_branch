#ifndef _LOGLIKE_THETA_FAST_H_
#define _LOGLIKE_THETA_FAST_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void 	loglike_theta_fast_initialize();
void 	loglike_theta_fast_finalize();
float 	loglike_theta_fast(float *x, int n, void *output );

float 	loglike_(float *x, int n, void *output, int *info );


#endif
