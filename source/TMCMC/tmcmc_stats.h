
#ifndef _TMCMC_STATS_H_
#define _TMCMC_STATS_H_

void calculate_statistics(float flc[], unsigned int n, int nselections, int gen, unsigned int sel[]);


typedef struct fparam_s {
	float *fj;
	int     fn;
	float  pj;
	float  tol;
} fparam_t;



#endif
