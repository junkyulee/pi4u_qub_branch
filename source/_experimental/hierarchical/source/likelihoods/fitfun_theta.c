#include "../TMCMC/fitfun.h"
#include "loglike_theta.h"


void fitfun_initialize() {

}



void fitfun_finalize() {

}




float fitfun(float *x, int n, void *output, int *info) {

	return loglike_theta( x, n, output, info );

}
