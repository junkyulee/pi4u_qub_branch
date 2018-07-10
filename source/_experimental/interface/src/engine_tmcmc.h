/*
 *  engine_tmcmc.h
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
//#include <mpi.h>
//#include <torc.h>

#include "gsl_headers.h"

typedef struct data_s {
	int	Nth;		// = PROBDIM
	int	MaxStages;	// = MAXGENS
	int	PopSize;	// = DATANUM

	float	*lowerbound;	//[PROBDIM];
	float	*upperbound;	//[PROBDIM];

	float lb, ub;		//generic lower and upper bound

	float	TolCOV;
	float	bbeta;
	int	seed;

	struct optim_options {
		int	MaxIter;
		float	Tol;
		int	Display;
	} options;

	int	iplot;
	int	idump;

	int	*Num;		//[MAXGENS];
	int	LastNum;

} data_t;

typedef struct runinfo_s {
	int 	Gen;
	float	*CoefVar;		//[MAXGENS];
	float	*p;			//[MAXGENS];		// cluster-wide
	int	*currentuniques;	//[MAXGENS];
	float	*logselection;		//[MAXGENS];
	float	*acceptance;		//[MAXGENS];
	float	**SS;			//[PROBDIM][PROBDIM];	// cluster-wide
	float	**meantheta; 		//[MAXGENS][PROBDIM]
} runinfo_t;


/*** DATABASES ***/
typedef struct cgdbp_s {
	float *point; //[PROBDIM];
	float F;
	int counter;	// not used (?)
	int nsel;	// for selection of leaders only
	int queue;	// for submission of leaders only
} cgdbp_t;

typedef struct cgdb_s {
	cgdbp_t *entry; //[MAX_DB_ENTRIES];
	int entries;
#if 0
	pthread_mutex_t m;
#endif
} cgdb_t;


typedef struct tmcmc_data_s {
	data_t data;
	runinfo_t runinfo;
	cgdb_t curgen_db;
} tmcmc_data_t;

extern tmcmc_data_t tmcmc_data;
//extern data_t data;
//extern runinfo_t runinfo;
//extern cgdb_t curgen_db;

void update_curgen_db(cgdb_t *curgen_db, float point[], float F, int dim);
void init_curgen_db(cgdb_t *curgen_db, int popsize);
void dump_curgen_db(cgdb_t *curgen_db, int Gen, int dim, int tinfo[4]);

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
void shuffle(int *perm, int N);
int mvnrnd(float *mean, float *var, float *res, int n);

/*** STATISTICS ***/
void calculate_statistics(tmcmc_data_t *tmcmc_data, float flc[], int n, int nselections, int gen, unsigned int sel[]);

/*** PROBLEM FUNCTIONS ***/
float likelihood(float *x, int N);
float posterior(float *theta, int n, float LH);
float priorpdf(float *theta, int n, float *lowerbound, float *upperbound);

/*** AUX ***/
void inc_nfc(void);
void get_nfc_task(int *);
int get_nfc(void);
void reset_nfc_task(void);
void reset_nfc(void);
int get_tfc(void);


int tmcmc_initialize(char *fitfun_name);
void tmcmc_finalize(void);
void tmcmc(float *res, int tmcmc_info[4], int Nth, int MaxStages, int PopSize, float *lb, float *ub);
