#ifndef _FITFUN_H_
#define _FITFUN_H_

void fitfun_initialize(char *name);
void fitfun_finalize();
float fitfun(float *x, int N, void *output, int *info);

#endif
