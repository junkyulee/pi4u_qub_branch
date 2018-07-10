/*
 *  problem.h
 *  Pi4U
 *
 *  Created by Lina Kulakova on 1/1/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */


#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_float.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <mutex>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "fitfun.h"
#include "random.h"


// dimensionality of the problem (including errors)
int PROBDIM=2;	// user input

// probabilistic ABC; if not defined, use standard ABC
//#define PROBABILISTIC

int gsl_vector_is_positive(const gsl_vector * v);
void gsl_vector_fprintf_in_line(FILE * fp, const gsl_vector * v, const char * format);

// generate only positive parameters
int enforce_positivity[] = {0, 0, 0};		// user input
float lower_bounds[] = {-6.0, -6.0, -6.0};	// user input
float upper_bounds[] = {+6.0, +6.0, +6.0};	// user input
float normalization[] = {1, 1, 1};		// user input

float eval_count = 0;

void generate_sample_by_prior(gsl_vector * sample)
{
	do
	{
		// GSL is thread-safe
		// generate each coordinate according to the prior distribution
		for (int i = 0; i < PROBDIM; i++) {
			gsl_vector_set(sample, i, gsl_ran_flat(rng.r[torc_i_worker_id()], lower_bounds[i]*normalization[i], upper_bounds[i]*normalization[i]));
		}
	}
	while(!gsl_vector_is_positive(sample));
}

float evaluate_prior(const gsl_vector * sample)
{
	float res = 1.0;

	for (int i = 0; i < PROBDIM; i++) {
		res *= gsl_ran_flat_pdf(gsl_vector_get(sample, i), lower_bounds[i]*normalization[i], upper_bounds[i]*normalization[i]);
	}

	return res;
}

float evaluate_logprior(const gsl_vector * sample, int * isinf)
{
	float res = 0.0;

	for (int i = 0; i < PROBDIM; i++) {
		res += gsl_ran_flat_logpdf(gsl_vector_get(sample, i), lower_bounds[i]*normalization[i], upper_bounds[i]*normalization[i], isinf);
		if (*isinf) return 0.0;	// peh: check this
	}

	return res;
}

void run_problem(float * sample_v, float * discrepancy, int *info)
{
	int n = PROBDIM;

	*discrepancy = fitfun(sample_v, n, NULL, info);
}

void normalize_sample(gsl_vector * sample)
{
	for(int i=0; i<sample->size; ++i)
		gsl_vector_set(sample, i, gsl_vector_get(sample, i)/normalization[i]);
}
