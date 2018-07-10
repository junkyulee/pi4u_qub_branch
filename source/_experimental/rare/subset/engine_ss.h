/*
 *  engine_ss.h
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <torc.h>

/*
#ifndef MY_GETTIME
#define MY_GETTIME
#include <sys/time.h>
static float my_gettime()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return (float)t.tv_sec + (float)t.tv_usec*1.0E-6;
}
#endif
*/

#include "gsl_headers.h"

typedef struct data_s {
	int	N_init;		//  = 1100*10;
	int	N_seeds;	// = 100*10;
	int	N_steps;	// = 10;

	int	Nth;	// = 2;

	float	*lowerbound;	//[PROBDIM] = {-6.0, -6.0};
	float	*upperbound;	//[PROBDIM] = {+6.0, +6.0};
	float	lb, ub;		// generic lower and upper bound

	int	NTHRESHOLDS;	//	= 10;
	float	*threshold;	// [NTHRESHOLDS] = {-64.0, -32.0, -16.0, -8.0, -4.0, -2.0, -1.0, -1.0/2.0, -1.0/4.0, -1.0/8.0};

	int	logval;

	float	sigma;
	int	seed;
	int	iplot;


	float	FACTOR;		// for adaptive subset
	int	MAXTHRESHOLDS;	//

} data_t;


//typedef struct runinfo_s {
//} runinfo_t;

extern data_t data;
//extern runinfo_t runinfo;


/*** DATABASES ***/

void db_init();
void add_seed_task(float s[], float *pfs);
void add_seed(float s[], float *pfs);
void add_sample_task(float s[], float *pfs);
void add_sample(float s[], float *pfs);
void permute_seeds(int N);
int compar_desc(const void* p1, const void* p2);
int compar_asc(const void* p1, const void* p2);
void sort_seeds(int N, int desc);
void dump_samples(int step);
void dump_seeds(int step);

void set_nsamples(int val);
void set_nseeds(int val);
int get_nsamples();
int get_nseeds();
float *get_sample(int i);
float *get_sample_f_addr(int i);
float get_sample_f(int i);
float *get_seed(int i);
float get_seed_f(int i);
float *get_seed_f_addr(int i);

/*** UTILS ***/
float compute_sum(float *x, int n);
float compute_mean(float *x, int n);
float compute_std(float *x, int n, float mean);
float compute_min(float *x, int n);
int compute_min_idx_i(int *v, int n);
float compute_max(float *x, int n);
void print_matrix(char *name, float *x, int n);
void print_matrix_i(char *name, int *x, int n);
void print_matrix_2d(char *name, float **x, int n1, int n2);

/*** RNG ***/
void gsl_rand_init(int seed);
float normalrand(float mu, float var);
float uniformrand(float a, float b);
void multinomialrand(size_t K, unsigned int N, float q[], unsigned int nn[]);
void shuffle(int *data, int N);
int mvnrnd(float *mean, float *var, float *res, int n);

/*** AUX ***/
void inc_nfc();
void get_nfc_task(int *);
int get_nfc();
void reset_nfc_task();
void reset_nfc();
int get_tfc();
