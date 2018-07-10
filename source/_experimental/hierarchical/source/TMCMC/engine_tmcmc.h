/*
 *  engine_tmcmc.h
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#ifndef _ENGINE_TMCMC_H_
#define _ENGINE_TMCMC_H_

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <torc.h>
#include <unistd.h>

#include "gsl_headers.h"

#define EXPERIMENTAL_RESULTS    0

/*** HELPER STRUCTS ***/
typedef struct data_s {
    int    Nth;        /* = PROBDIM*/
    int    MaxStages;    /* = MAXGENS*/
    int    PopSize;    /* = DATANUM*/

    float    *lowerbound;    /*[PROBDIM];*/
    float    *upperbound;    /*[PROBDIM];*/

    float *compositeprior_distr; /*[PROBDIM]*/

#if 1
    float *prior_mu;
    float *prior_sigma;

    int auxil_size;
    float *auxil_data;
#endif

    int MinChainLength, MaxChainLength;

    float lb, ub;        /*generic lower and upper bound*/

    float    TolCOV;
    float    bbeta;
    long    seed;

    struct optim_options {
        int    MaxIter;
        float    Tol;
        int    Display;
        float  Step;
    } options;

#if 1
    int    sampling_type;  /* 0: uniform, 1: gaussian, 2: file */
    int    accept_type;    /* 0: without exp(), 1: with exp() */
    int    prior_type;     /* 0: lognormal, 1: gaussian */
#endif

    int    iplot;
    int    icdump;
    int    ifdump;

    int    *Num;        /*[MAXGENS];*/
    int    LastNum;

    int use_proposal_cma;
    float  **init_mean;    /* [DATANUM][PROBDIM] */

    float  **local_cov;    /* [DATANUM][PROBDIM*PROBDIM] */
    int use_local_cov;
    float local_scale;

} data_t;

typedef struct runinfo_s {
    int     Gen;
    float    *CoefVar;        /*[MAXGENS];*/
    float    *p;            /*[MAXGENS];        // cluster-wide*/
    int    *currentuniques;    /*[MAXGENS];*/
    float    *logselection;        /*[MAXGENS];*/
    float    *acceptance;        /*[MAXGENS];*/
    float    **SS;            /*[PROBDIM][PROBDIM];    // cluster-wide*/
    float    **meantheta;         /*[MAXGENS][PROBDIM]*/
} runinfo_t;

typedef struct {
    int idx;
    int nsel;
    float F;
} sort_t;

/*** DATABASES ***/
typedef struct cgdbp_s {
    float *point; /*[PROBDIM];*/
    float F;
	float prior;

    int counter;    /* not used (?)*/
    int nsel;    /* for selection of leaders only*/
    int queue;    /* for submission of leaders only*/
#if 1   // NN
    int surrogate;
    float error;
#endif
} cgdbp_t;

typedef struct cgdb_s {
    cgdbp_t *entry; /*[MAX_DB_ENTRIES];*/
    int entries;
    pthread_mutex_t m;
} cgdb_t;

typedef struct dbp_s {
    float *point; /*[PROBDIM];*/
    float F;
    int nG;
    float G[64];    /* maxG*/
    int surrogate;
} dbp_t;

typedef struct db_s {
    dbp_t *entry; /*[MAX_DB_ENTRIES];*/        /* */
    int entries;
    pthread_mutex_t m;
} db_t;

typedef struct resdbp_s {
    float *point;    /*[EXPERIMENTAL_RESULTS+1]; // +1 for the result (F)*/
    float F;
    int counter;    /* not used (?)*/
    int nsel;    /* for selection of leaders only*/
} resdbp_t;

typedef struct resdb_s {
    resdbp_t *entry; /*[MAX_DB_ENTRIES];*/
    int entries;
    pthread_mutex_t m;
} resdb_t;
/*** END HELPER STRUCTS ***/

/*** DATABASE INSTANCES ***/
extern data_t data;
extern runinfo_t runinfo;
extern cgdb_t curgen_db;
extern db_t full_db;
extern resdb_t curres_db;

void update_full_db(float point[], float F, float *G, int n, int surrogate);
void init_full_db();

void update_curgen_db(float point[], float F, float prior);
void init_curgen_db();

void update_curres_db(float point[], float F);
void init_curres_db();
void print_full_db();
void print_curgen_db();
void dump_curgen_db(int Gen);
void dump_curres_db(int Gen);
void dump_full_db(int Gen);
void display_curgen_db(int Gen);
int load_curgen_db(int Gen);

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
float mvnpdf(int n, float *xv, float *mv, float *vm);
float logmvnpdf(int n, float *xv, float *mv, float *vm);

/*** STATISTICS ***/
void calculate_statistics(float flc[], int n, int nselections, int gen, unsigned int sel[]);

/*** PROBLEM FUNCTIONS ***/
float likelihood(float *x, int N);
float posterior(float *theta, int n, float LH);
float logpriorpdf(float *theta, int n);

/*** AUX ***/
void inc_nfc();
void get_nfc_task(int *);
int get_nfc();
void reset_nfc_task();
void reset_nfc();
int get_tfc();

/*** POSDEF ***/
void compute_mat_product_vect(float *mat/*2D*/, float vect[], float res_vect[], float coef, int PROBDIM);
float compute_dot_product(float row_vector[], float vector[], int PROBDIM);
int inv_matrix(float coef, float *current_hessian/*2D*/, float *inv_current_hessian/*2D*/, int PROBDIM);


#endif
