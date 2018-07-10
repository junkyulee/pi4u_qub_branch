/*
 *  tmcmc_engine.h
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#ifndef _ENGINE_TMCMC_H_
#define _ENGINE_TMCMC_H_




#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>



#if defined(_USE_TORC_)
	#include <mpi.h>
	#ifdef __cplusplus
		extern "C"
		{
	#endif

	#include <torc.h>

	#ifdef __cplusplus
		}
	#endif

#else
	#include <pthread.h>
#endif



#define EXPERIMENTAL_RESULTS    0






/*** HELPER STRUCTS ***/
typedef struct data_s {
    int    Nth;        /* = PROBDIM*/
    int    MaxStages;    /* = MAXGENS*/
    int    PopSize;    /* = DATANUM*/

    float    *lowerbound;    /*[PROBDIM];*/
    float    *upperbound;    /*[PROBDIM];*/

    float *compositeprior_distr; /*[PROBDIM]*/

    float *prior_mu;
    float *prior_sigma;

    int auxil_size;
    float *auxil_data;

    int MinChainLength, MaxChainLength;

    float lb, ub;        /*generic lower and upper bound*/

    float    TolCOV;
    float    bbeta;
    long    seed;
    int    burn_in;

    struct optim_options {
        int    MaxIter;
        float    Tol;
        int    Display;
        float  Step;
    } options;

    int    prior_type;     /* 0: uniform, 1: gaussian, 3: composite */
    int    load_from_file;

    int    icdump;
    int    ifdump;

    int    *Num;        /*[MAXGENS];*/
    int    LastNum;

    int use_proposal_cma;
    float  **init_mean;    /* [DATANUM][PROBDIM] */

    float  **local_cov;    /* [DATANUM][PROBDIM*PROBDIM] */
    int use_local_cov;
    float local_scale;

    int stealing;
    int restart;
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
    int surrogate;
    float error;
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
extern data_t		data;
extern runinfo_t 	runinfo;
extern cgdb_t 		curgen_db;
extern db_t 		full_db;
extern resdb_t 		curres_db;






void read_data();
void data_init();
void save_runinfo();
int load_runinfo();
void check_for_exit();
void torc_update_full_db_task(float point[], float *pF, int *psurrogate);
void torc_update_full_db_task_p5(float point[], float *pF, float *G, int *pn, int *psurrogate);
void torc_update_full_db(float point[], float F, float *G, int n, int surrogate);
void torc_update_curgen_db_task(float point[], float *pF, float *pprior);
void torc_update_curgen_db(float point[], float F, float prior);
void torc_update_curres_db_task(float point[EXPERIMENTAL_RESULTS], float *pF);
void torc_update_curres_db(float point[EXPERIMENTAL_RESULTS], float F);

/*** TASK MANAGEMENT ***/
void taskfun(float /*const*/ *x, int *pN, float *res, int winfo[4]);
float F(float *TP, int *pn);    /* for PNDL */
void evaluate_F(float point[], float *Fval, int worker_id, int gen_id, int chain_id, int step_id, int ntasks);


void initchaintask(float in_tparam[], int *pdim, float *out_tparam, int winfo[4]);



static int in_rect(float *v1, float *v2, float *diam, float sc, int D);
void precompute_chain_covariances(const cgdbp_t* leader,float** init_mean, float** chain_cov, int newchains);
int compute_candidate(float candidate[], float chain_mean[], float var);
int compute_candidate_cov(float candidate[], float chain_mean[], float chain_cov[]);
void chaintask(float in_tparam[], int *pdim, int *pnsteps, float *out_tparam, int winfo[4],float *init_mean, float *chain_cov);
int compar_desc(const void* p1, const void* p2);
int prepare_newgen(int nchains, cgdbp_t *leaders);







#endif
