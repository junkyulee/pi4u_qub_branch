/*
 *  dram.h
 *  Pi4U
 *
 *  Copyright 2018 ETH Zurich. All rights reserved.
 *
 */

#ifndef _DRAM_H_
#define _DRAM_H_

#include "../priors/priors.h"



/* Variables related to the options of DRAM */
struct _params
{
	int 	Npar;    	/* problem dimensionality (number of parameters) */
    int 	Nsim;    	/* Length of the chain */

    float 	DRscale;    /* DR shrink factor, if zero, no DR */

	int 	AMinterv;	/* how often to adapt, if zero, no adaptation */	// 20 //= 500;
	float 	AMscale;	/* Scale for adapting the proposal = 2.4/sqrt(Npar) */
	float 	AMepsilon;	/* Small constant used in the AM covariance update */

	float 	*qcov;     	/* Proposal covariance matrix */

    int     seed;

	int 	printfreq;	/* Time interval for printing info on screen */
	int 	verbose;


	float 	*par0;		/* Initial parameter vector */
	float 	sigma2;		/* Initial/prior value for the Gaussian error variance */
	float 	n0;			/* Precision of sigma2 as imaginative observations. if n0<0, no sigma2 update */
	int 	n;			/* Nnumber of actual observations (for sigma2 update) */

	float 	*lbounds;	/* Lower bounds of parameter values */
	float 	*ubounds;	/* Upper bounds of parameter values */

	Density *prior;

	char 	filename[256];	/* Output file for the chain */

} par;



void dram_init();
void dram();
void dram_finalize();



#endif
