/*
 *  engine_tmcmc.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <assert.h>
#include <signal.h>
#include <stdio.h>
#include "engine_tmcmc.h"
#include "fitfun.h"

/*#define _STEALING_*/
/*#define VERBOSE 1*/
/*#define _RESTART_*/

#if !defined(_USE_TORC_)
static int torc_node_id() { return 0; }

#if defined(_USE_OPENMP_)
#include <omp.h>
static int torc_worker_id() { return omp_get_thread_num(); }
#else
static int torc_worker_id() { return 0; }
#endif

#include <sys/time.h>
static float torc_gettime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (float)t.tv_sec + (float)t.tv_usec*1.0E-6;
}
#endif

data_t data;
runinfo_t runinfo;
cgdb_t curgen_db;
db_t full_db;
resdb_t curres_db;

void read_data()
{
    int i, found;

    /* DEFAULT VALUES */
    data.Nth = 4; /* Default PROBDIM */
    data.MaxStages = 20; /* Default MAXGENS; */
    data.PopSize = 1024;    /* Default DATANUM */

    data.MinChainLength = 0;
    data.MaxChainLength = 1e6;

    data.lb = -6;    /* Default LB, same for all */
    data.ub = +6;    /* Default UB, same for all */

    data.lowerbound = (float *)malloc(data.Nth*sizeof(float));
    data.upperbound = (float *)malloc(data.Nth*sizeof(float));

    for (i = 0; i < data.Nth; i++) {
        data.lowerbound[i] = data.lb;
        data.upperbound[i] = data.ub;
    }

    /* new */
    data.prior_mu = (float *)malloc(data.Nth*sizeof(float));
    data.prior_sigma = (float *)malloc(data.Nth*data.Nth*sizeof(float));
    data.auxil_size = 0;
    data.auxil_data = NULL;

    for (int i = 0; i < data.Nth; i++) {
        data.prior_mu[i] = 0.0;
    }

    for (i = 0; i < data.Nth; i++) {
        int j;
        for (j = 0; j < data.Nth; j++) {
            if (i == j)
                data.prior_sigma[i*data.Nth+j] = 1.0;
            else
                data.prior_sigma[i*data.Nth+j] = 0.0;
        }
    }

    data.TolCOV = 1.0;    /* 0.25, 0.5 */
    data.bbeta = 0.2;
    data.seed = 280675;
    data.burn_in = 0;

    data.options.MaxIter = 100;    /**/
    data.options.Tol = 1e-6;
    data.options.Display = 0;
    data.options.Step = 1e-5;

    data.prior_type = 0;    /* uniform > gaussian > composite */
    data.load_from_file = 0;    /* load initial samples from file instead from prior */

    data.icdump = 1;    /* dump current dataset of accepted points */
    data.ifdump = 0;    /* dump complete dataset of points */

    data.Num = (int *)malloc(data.MaxStages*sizeof(int));
    for (i = 0; i < data.MaxStages; i++) {
        data.Num[i] = data.PopSize; /* default DATANUM */
    }
    data.LastNum = data.PopSize; /* DATANUM; */

    data.stealing = 0;
    data.restart = 0;

    /* USER-DEFINED VALUES */
    FILE *f = fopen("tmcmc.par", "r");
    if (f == NULL) {
        return;
    }

    char line[256];

    int line_no = 0;
    while (fgets(line, 256, f)!= NULL) {
        line_no++;
        if ((line[0] == '#')||(strlen(line)==0)) {
#if VERBOSE
            printf("ignoring line %d\n", line_no);
#endif
            continue;
        }

        if (strstr(line, "Nth")) {
            sscanf(line, "%*s %d", &data.Nth);
        }
        else if (strstr(line, "MaxStages")) {
            sscanf(line, "%*s %d", &data.MaxStages);
        }
        else if (strstr(line, "PopSize")) {
            sscanf(line, "%*s %d", &data.PopSize);
        }
        else if (strstr(line, "TolCOV")) {
            sscanf(line, "%*s %lf", &data.TolCOV);
        }
        else if (strstr(line, "bbeta")) {
            sscanf(line, "%*s %lf", &data.bbeta);
        }
        else if (strstr(line, "seed")) {
            sscanf(line, "%*s %ld", &data.seed);
        }
        else if (strstr(line, "burn_in")) {
            sscanf(line, "%*s %d", &data.burn_in);
        }
        else if (strstr(line, "opt.MaxIter")) {
            sscanf(line, "%*s %d", &data.options.MaxIter);
        }
        else if (strstr(line, "opt.Tol")) {
            sscanf(line, "%*s %lf", &data.options.Tol);
        }
        else if (strstr(line, "opt.Display")) {
            sscanf(line, "%*s %d", &data.options.Display);
        }
        else if (strstr(line, "opt.Step")) {
            sscanf(line, "%*s %lf", &data.options.Step);
            printf("setting step = %f\n", data.options.Step);
        }
        else if (strstr(line, "prior_type")) {
            sscanf(line, "%*s %d", &data.prior_type);
        }
        else if (strstr(line, "load_from_file")) {
            sscanf(line, "%*s %d", &data.load_from_file);
        }
        else if (strstr(line, "icdump")) {
            sscanf(line, "%*s %d", &data.icdump);
        }
        else if (strstr(line, "ifdump")) {
            sscanf(line, "%*s %d", &data.ifdump);
        }
        else if (strstr(line, "Bdef")) {
            sscanf(line, "%*s %lf %lf", &data.lb, &data.ub);
        }
        else if (strstr(line, "MinChainLength")) {
            sscanf(line, "%*s %d", &data.MinChainLength);
        }
        else if (strstr(line, "MaxChainLength")) {
            sscanf(line, "%*s %d", &data.MaxChainLength);
        }
        else if (strstr(line, "use_local_cov")) {
            sscanf(line, "%*s %d", &data.use_local_cov);
        }
        else if (strstr(line, "stealing")) {
            sscanf(line, "%*s %d", &data.stealing);
        }
        else if (strstr(line, "restart")) {
            sscanf(line, "%*s %d", &data.restart);
        }
    }

    rewind(f);
    line_no = 0;

    free(data.lowerbound);
    free(data.upperbound);
    data.lowerbound = (float *)malloc(data.Nth*sizeof(float));
    data.upperbound = (float *)malloc(data.Nth*sizeof(float));

    /* Read the lower and upper bounds for each parameter - used by uniform and truncated gaussian priors */
    for (i = 0; i < data.Nth; i++) {
        found = 0;
        while (fgets(line, 256, f)!= NULL) {
            line_no++;

            if ((line[0] == '#')||(strlen(line)==0)) continue;

            char bound[16];
            sprintf(bound, "B%d", i);
            if (strstr(line, bound) != NULL) {
                sscanf(line, "%*s %lf %lf", &data.lowerbound[i], &data.upperbound[i]);
                found = 1;
                break;
            }
        }
        if (!found) {
            data.lowerbound[i] = data.lb;    /* Bdef value or Default LB */
            data.upperbound[i] = data.ub;    /* Bdef value of Default UB */
        }
        rewind(f);
        line_no = 0;
    }

    if(data.prior_type == 1) /* gaussian */
    {
        /* new, parse prior_mu */
        rewind(f);
        line_no = 0;

        free(data.prior_mu);
        data.prior_mu = (float *)malloc(data.Nth*sizeof(float));

        found = 0;
        while (fgets(line, 256, f)!= NULL) {
            line_no++;
            if ((line[0] == '#')||(strlen(line)==0)) continue;

            if (strstr(line, "prior_mu") != NULL) {
                char *tok = strtok(line, " ;,\t");
                if (tok == NULL) break;
                int i = 0;
                tok = strtok(NULL, " ;,\t");
                while (tok != NULL) {
                    data.prior_mu[i] = atof(tok);
                    i++;
                    tok = strtok(NULL, " ;,\t");
                }
                found = 1;
                break;
            }
        }

        if (!found) {
            for (i = 0; i < data.Nth; i++) {
                data.prior_mu[i] = 0.0;        /* Mudef value of Default Mean */
            }
        }

        /* new, parse prior_sigma */
        rewind(f);
        line_no = 0;

        free(data.prior_sigma);
        data.prior_sigma = (float *)malloc(data.Nth*data.Nth*sizeof(float));

        found = 0;
        while (fgets(line, 256, f)!= NULL) {
            line_no++;
            if ((line[0] == '#')||(strlen(line)==0)) continue;

            if (strstr(line, "prior_sigma") != NULL) {
                char *tok = strtok(line, " ;,\t");
                if (tok == NULL) break;
                int i = 0;
                tok = strtok(NULL, " ;,\t");
                while (tok != NULL) {
                    data.prior_sigma[i] = atof(tok);
                    i++;
                    tok = strtok(NULL, " ;,\t");
                }
                found = 1;
                break;
            }
        }

        if (!found) {
            for (i = 0; i < data.Nth; i++) {
                int j;
                for (j = 0; j < data.Nth; j++) {
                    if (i == j)
                        data.prior_sigma[i*data.Nth+j] = 1.0;    /* Sigmadef value of Default Sigma */
                    else
                        data.prior_sigma[i*data.Nth+j] = 0.0;
                }
            }
        }
    }

    if(data.prior_type == 3) /* composite */
    {
        rewind(f);
        line_no = 0;

        data.compositeprior_distr = (float *)malloc(data.Nth*sizeof(float));

        free(data.prior_mu);
        free(data.prior_sigma);
        data.prior_mu = (float *)malloc(data.Nth*sizeof(float));
        data.prior_sigma = (float *)malloc(data.Nth*data.Nth*sizeof(float));

        for (i = 0; i < data.Nth; i++) {
            found = 0;
            while (fgets(line, 256, f)!= NULL) {
                line_no++;

                if ((line[0] == '#')||(strlen(line)==0)) continue;

                char bound[16];
                sprintf(bound, "C%d", i);
                if (strstr(line, bound) != NULL) {
                    float type, param0, param1;
                    sscanf(line, "%*s %lf %lf %lf", &type, &param0, &param1);

                    printf("type = %lf, x0 = %lf, x1 = %lf\n", type, param0, param1);

                    data.compositeprior_distr[i] = type;
                    if (type == 0) {
                        data.lowerbound[i] = param0;
                        data.upperbound[i] = param1;
                    }
                    if ((type == 1) || (type == 2)) {
                        data.prior_mu[i] = param0;
                        data.prior_sigma[i] = param1;
                    }

                    found = 1;
                    break;
                }
            }
            if (!found) {
                data.lowerbound[i] = data.lb;    /* Bdef value or Default LB */
                data.upperbound[i] = data.ub;    /* Bdef value of Default UB */
                data.compositeprior_distr[i] = 0;
            }
            rewind(f);
            line_no = 0;
        }
    }

    /* new, parse auxil_size and auxil_data */
    rewind(f);
    line_no = 0;

    found = 0;
    while (fgets(line, 256, f)!= NULL) {
        line_no++;
        if ((line[0] == '#')||(strlen(line)==0)) continue;

        if (strstr(line, "auxil_size") != NULL) {
            sscanf(line, "%*s %d", &data.auxil_size);
            found = 1;
            break;
        }
    }

    if (data.auxil_size > 0)
    {
        rewind(f);
        line_no = 0;

        data.auxil_data = (float *)malloc(data.auxil_size*sizeof(float));

        found = 0;
        while (fgets(line, 256, f)!= NULL) {
            line_no++;
            if ((line[0] == '#')||(strlen(line)==0)) continue;

            if (strstr(line, "auxil_data") != NULL) {
                char *tok = strtok(line, " ;,\t");
                if (tok == NULL) break;
                int i = 0;
                tok = strtok(NULL, " ;,\t");
                while (tok != NULL) {
                    data.auxil_data[i] = atof(tok);
                    i++;
                    tok = strtok(NULL, " ;,\t");
                }
                found = 1;
                break;
            }
        }
    }

    fclose(f);


#if 0
    print_matrix((char *)"prior_mu", data.prior_mu, data.Nth);
    print_matrix((char *)"prior_sigma", data.prior_sigma, data.Nth*data.Nth);
    print_matrix((char *)"auxil_data", data.auxil_data, data.auxil_size);
#endif

    free(data.Num);
    data.Num = (int *)malloc(data.MaxStages*sizeof(int));
    for (i = 0; i < data.MaxStages; i++) {
        data.Num[i] = data.PopSize;
    }
    data.LastNum = data.PopSize;

    float *LCmem = (float *)calloc(1, data.PopSize*data.Nth*data.Nth*sizeof(float));
    data.local_cov = (float **)malloc(data.PopSize*sizeof(float *));
    int pos;
    for (pos = 0; pos < data.PopSize; ++pos)
    {
        data.local_cov[pos] = LCmem + pos*data.Nth*data.Nth;
        for (i=0; i<data.Nth; ++i)
            data.local_cov[pos][i*data.Nth+i] = 1;
    }
}

void data_init()
{
    int i;

    /* DATA: user's input parameters */
    read_data();

#if 1
    init_curgen_db();
    init_curres_db();
    init_full_db();
#endif

    /* RUNINFO: running state */
    runinfo.CoefVar = (float *)calloc(1, (data.MaxStages+1)*sizeof(float));
    runinfo.p = (float *)calloc(1, (data.MaxStages+1)*sizeof(float));
    runinfo.currentuniques = (int *)calloc(1, data.MaxStages*sizeof(int));
    runinfo.logselection = (float *)calloc(1, data.MaxStages*sizeof(float));
    runinfo.acceptance = (float *)calloc(1, data.MaxStages*sizeof(float));

    float *SSmem = (float *)calloc(1, data.Nth*data.Nth*sizeof(float));
    runinfo.SS = (float **)malloc(data.Nth*sizeof(float *));
    for (i = 0; i < data.Nth; i++) {
        runinfo.SS[i] = SSmem + i*data.Nth; /*&SSmem[i*data.Nth];*/
    }

    runinfo.meantheta = (float **)calloc(1, (data.MaxStages+1)*sizeof(float *));
    for (i = 0; i < data.MaxStages+1; i++) {
        runinfo.meantheta[i] = (float *)calloc(1, data.Nth*sizeof(float));
    }

    runinfo.Gen = 0;
    runinfo.CoefVar[0] = 10;

    printf("runinfo = %p\n", &runinfo);
    printf("runinfo.p = %p\n", runinfo.p);
    printf("runinfo.SS = %p\n", runinfo.SS);
}

void save_runinfo()
{
    int i, j;

    /* allocate and initialize runinfo */
    FILE *f = fopen("runinfo.txt", "w");

    fprintf(f, "Gen=\n");
    fprintf(f, "%d\n", runinfo.Gen);

    fprintf(f, "CoefVar[%d]=\n", data.MaxStages);
    for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.CoefVar[i]);

    fprintf(f, "p[%d]=\n", data.MaxStages);
    for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.p[i]);

    fprintf(f, "currentuniques[%d]=\n", data.MaxStages);
    for (i = 0; i < data.MaxStages; i++) fprintf(f, "%d\n", runinfo.currentuniques[i]);

    fprintf(f, "logselection[%d]=\n", data.MaxStages);
    for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.logselection[i]);

    fprintf(f, "acceptance[%d]=\n", data.MaxStages);
    for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.acceptance[i]);

    fprintf(f, "SS[%d][%d]=\n", data.Nth, data.Nth);
    for (i = 0; i < data.Nth; i++)
        for (j = 0; j < data.Nth; j++)
            fprintf(f, "%.16lf\n", runinfo.SS[i][j]);

    fprintf(f, "meantheta[%d][%d]\n", data.MaxStages, data.Nth);
    for (i = 0; i < data.MaxStages; i++)
        for (j = 0; j < data.Nth; j++)
            fprintf(f, "%.16lf\n", runinfo.meantheta[i][j]);

    fclose(f);
}

int load_runinfo()
{
    int i, j;
    char header[256];

    /* allocate and initialize runinfo */
    FILE *f = fopen("runinfo.txt", "r");
    if (f == NULL) return 1;

    fscanf(f, "%s", header);
    fscanf(f, "%d", &runinfo.Gen);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.CoefVar[i]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.p[i]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++) fscanf(f, "%d\n", &runinfo.currentuniques[i]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.logselection[i]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.acceptance[i]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.Nth; i++)
        for (j = 0; j < data.Nth; j++)
            fscanf(f, "%lf\n", &runinfo.SS[i][j]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++)
        for (j = 0; j < data.Nth; j++)
            fscanf(f, "%lf\n", &runinfo.meantheta[i][j]);

    fclose(f);

    return 0;
}

void check_for_exit()
{
    int val, exitgen = -1;
    char *s;

    s = (char *) getenv("EXIT_GEN");
    if (s != 0 && sscanf(s, "%d", &val) == 1 && val >= 0)
        exitgen = val;

    if (exitgen == runinfo.Gen) {
        printf("Read Exit Envrironment Variable!!!\n");
#if defined(_USE_TORC_)
        torc_finalize();
#endif
        exit(1);
    }

    FILE *fp;
    fp = fopen("exit.txt", "r");
    if (fp != NULL) {
        printf("Found Exit File!!!\n");
        unlink("exit.txt");
#if defined(_USE_TORC_)
        torc_finalize();
#endif
        exit(1);
    }
}

/* TORC-BASED DATA MANAGEMENT */
void torc_update_full_db_task(float point[], float *pF, int *psurrogate)
{
    float F = *pF;
    int surrogate = *psurrogate;
    float *G = NULL;
    int n = 0;

    update_full_db(point, F, G, n, surrogate);
}

void torc_update_full_db_task_p5(float point[], float *pF, float *G, int *pn, int *psurrogate)
{
    float F = *pF;
    int n = *pn;
    int surrogate = *psurrogate;

    update_full_db(point, F, G, n, surrogate);
}

void torc_update_full_db(float point[], float F, float *G, int n, int surrogate)
{
    if (torc_node_id() ==0) {
        update_full_db(point, F, G, n, surrogate);
        return;
    }

#if defined(_USE_TORC_)
    if (n == 0)
        torc_create_direct(0, (void (*)())torc_update_full_db_task, 3,        /* message to the database manager (separate process?) or direct execution by server thread */
                data.Nth, MPI_DOUBLE, CALL_BY_VAL,
                1, MPI_DOUBLE, CALL_BY_COP,
                1, MPI_INT, CALL_BY_COP,
                point, &F, &surrogate);
    else
        torc_create_direct(0, (void (*)())torc_update_full_db_task_p5, 5,
                data.Nth, MPI_DOUBLE, CALL_BY_VAL,
                1, MPI_DOUBLE, CALL_BY_COP,
                n, MPI_DOUBLE, CALL_BY_COP,    /* xxx: for CALL_BY_VAL: in the full-version of the library, with n=1 we had segv */
                1, MPI_INT, CALL_BY_COP,
                1, MPI_INT, CALL_BY_COP,
                point, &F, G, &n, &surrogate);

    torc_waitall3();
#endif
}

void torc_update_curgen_db_task(float point[], float *pF, float *pprior)
{
    float F = *pF;
	float prior = *pprior;

    update_curgen_db(point, F, prior);
}

void torc_update_curgen_db(float point[], float F, float prior)
{
    int me = torc_node_id();

    if (me == 0) {
        update_curgen_db(point, F,prior);
        return;
    }

#if defined(_USE_TORC_)
	torc_create_direct(0, (void (*)())torc_update_curgen_db_task, 3,            /* message to the database manager (separate process?) or direct execution by server thread */
		data.Nth, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_DOUBLE, CALL_BY_COP,
		point, &F,&prior);

    torc_waitall3();    /* wait without releasing the worker */
#endif
}

void torc_update_curres_db_task(float point[EXPERIMENTAL_RESULTS], float *pF)
{
    float F = *pF;

    update_curres_db(point, F);
}

void torc_update_curres_db(float point[EXPERIMENTAL_RESULTS], float F)
{
    int me = torc_node_id();

    if (me ==0) {
        update_curres_db(point, F);
        return;
    }

#if defined(_USE_TORC_)
    torc_create_direct(0, (void (*)())torc_update_curres_db_task, 2,        /* message to the database manager (separate process?) or direct execution by server thread */
            EXPERIMENTAL_RESULTS, MPI_DOUBLE, CALL_BY_COP,
            1, MPI_DOUBLE, CALL_BY_COP,
            point, &F);
    torc_waitall3();
#endif
}

/*** TASK MANAGEMENT ***/
void taskfun(float /*const*/ *x, int *pN, float *res, int winfo[4])
{
    float f;
    int N = *pN;

    inc_nfc();    /* increment function call counter*/

    f = fitfun(x, N, (void *)NULL, winfo);
#if (EXPERIMENTAL_RESULTS > 0)    /* peh: decide about this (results should be passed as argument to fitfun) */
    int i;
    float results[EXPERIMENTAL_RESULTS];
    for (i = 0; i < EXPERIMENTAL_RESULTS; i++) {
        if (i < data.Nth)
            results[i] = x[i];
        else
            results[i] = 0.0;
    }
    torc_update_curres_db(results, f);
#endif

    *res = f;
    return;
}

float F(float *TP, int *pn)    /* for PNDL */
{
    float gres;

    taskfun(TP, pn, &gres, NULL);

    return gres;
}


void evaluate_F(float point[], float *Fval, int worker_id, int gen_id, int chain_id, int step_id, int ntasks)
{
    float F;
    int winfo[4];
    int dim = data.Nth;

    winfo[0] = gen_id;
    winfo[1] = chain_id;
    winfo[2] = step_id;
    winfo[3] = 0;

#if VERBOSE
    printf("running on worker %d\n", worker_id);
#endif
    taskfun(point, &dim, &F, winfo);

    *Fval = F;
}

void initchaintask(float in_tparam[], int *pdim, float *out_tparam, int winfo[4])
{
    int i;
    int gen_id, chain_id;
    gen_id = winfo[0];
    chain_id = winfo[1];

    long me = torc_worker_id();
    float point[data.Nth], fpoint;

    for (i = 0; i < data.Nth; i++)
        point[i] = in_tparam[i];

    evaluate_F(point, &fpoint, me, gen_id, chain_id, 0, 1);

    float logprior = logpriorpdf(point, data.Nth);

    /* update current db entry */
    torc_update_curgen_db(point, fpoint, logprior);
    if (data.ifdump) torc_update_full_db(point, fpoint, NULL, 0, 0);
    *out_tparam = fpoint;    /* currently not required, the result is already in the db*/

    return;
}

static int in_rect(float *v1, float *v2, float *diam, float sc, int D)
{
    int d;
    for (d = 0; d < D; ++d) {
        if (fabs(v1[d]-v2[d]) > sc*diam[d]) return 0;
    }
    return 1;
}

void precompute_chain_covariances(const cgdbp_t* leader,
        float** init_mean, float** chain_cov, int newchains)
{
    printf("Precomputing covariances for the current generation...\n");

    int i, j, k, d, pos, ind;

    int D = data.Nth;
    int N = curgen_db.entries;

    float my_time = clock();

    // allocate space
    int* nn_ind = (int*)malloc(newchains*N*sizeof(int));
    int* nn_count = (int*)malloc(newchains*sizeof(int));
    float* diam = (float*)malloc(D*sizeof(float));
    float* chain_mean = (float*)malloc(D*sizeof(float));
    gsl_matrix* work = gsl_matrix_alloc(D, D);

    // find diameters
    for (d = 0; d < D; ++d) {
        float d_min = +1e6;
        float d_max = -1e6;
        for (pos = 0; pos < N; ++pos) {
            float s = curgen_db.entry[pos].point[d];
            if (d_min > s) d_min = s;
            if (d_max < s) d_max = s;
        }
        diam[d] = d_max-d_min;
        printf("Diameter %d: %.6lf\n", d, diam[d]);
    }

    int status = 0;
    float scale, ds = 0.05;
    for (scale = 0.1; scale <= 1.0; scale += ds) {
        // find neighbors in a rectangle - O(N^2)
        for (pos = 0; pos < newchains; ++pos) {
            nn_count[pos] = 0;
            float* curr = leader[pos].point;
            for (i = 0; i < N; ++i) {
                float* s = curgen_db.entry[i].point;
                if (in_rect(curr, s, diam, scale, D)) {
                    nn_ind[pos*N+nn_count[pos]] = i;
                    nn_count[pos]++;
                }
            }
        }

        // compute the covariances
        for (pos = 0; pos < newchains; ++pos) {
            for (d = 0; d < D; ++d) {
                chain_mean[d] = 0;
                for (k = 0; k < nn_count[pos]; ++k) {
                    ind = nn_ind[pos*N+k];
                    chain_mean[d] += curgen_db.entry[ind].point[d];
                }
                chain_mean[d] /= nn_count[pos];
            }

            for (i = 0; i < D; i++)
                for (j = 0; j < D; j++) {
                    float s = 0;
                    for (k = 0; k < nn_count[pos]; k++) {
                        ind = nn_ind[pos*N+k];
                        s += (curgen_db.entry[ind].point[i]-chain_mean[i]) *
                            (curgen_db.entry[ind].point[j]-chain_mean[j]);
                    }
                    chain_cov[pos][i*D+j] = chain_cov[pos][j*D+i] = s/nn_count[pos];
                }

            // check if the matrix is positive definite
            for (i = 0; i < D; ++i)
                for (j = 0; j < D; ++j) {
                    float s = chain_cov[pos][i*D+j];
                    gsl_matrix_set(work, i, j, s);
                }
            gsl_set_error_handler_off();
            status = gsl_linalg_cholesky_decomp(work);
            if (status == GSL_SUCCESS) break;
        }
    }

    if (status != GSL_SUCCESS) {
        for (i = 0; i < D; i++)
            for (j = 0; j < D; j++)
                chain_cov[pos][i*D+j] = data.bbeta*runinfo.SS[i][j];
    }

    free(nn_ind);
    free(nn_count);
    free(diam);
    free(chain_mean);
    gsl_matrix_free(work);

#if 0
    for (pos=0; pos<5; ++pos)
    {
        printf("Chain %d of %d: ", pos, newchains);
        print_matrix((char *)"chain_covariance", chain_cov[pos], D*D);
    }
#endif

    my_time = (clock() - my_time) / CLOCKS_PER_SEC;
    printf("Covariance computation time: %.2lf sec\n", my_time);
}

int compute_candidate(float candidate[], float chain_mean[], float var)
{
    int i, j;
    float bSS[data.Nth*data.Nth];

    for (i = 0; i < data.Nth; i++)
        for (j = 0; j < data.Nth; j++)
            bSS[i*data.Nth+j]= data.bbeta*runinfo.SS[i][j];

    //retry:
    mvnrnd(chain_mean, (float *)bSS, candidate, data.Nth);
    for (i = 0; i < data.Nth; i++) {
        if (isnan(candidate[i])) {
            printf("!!!!  isnan in candidate point!\n");
            exit(1);
            break;
        }
        if ((candidate[i] < data.lowerbound[i])||(candidate[i] > data.upperbound[i])) break;
    }
    //    if (i < data.Nth) goto retry;
    if (i < data.Nth) return -1;
    return 0;    // all good
}

int compute_candidate_cov(float candidate[], float chain_mean[], float chain_cov[])
{
    int i;

    mvnrnd(chain_mean, (float *)chain_cov, candidate, data.Nth);
    for (i = 0; i < data.Nth; i++) {
        if (isnan(candidate[i])) {
            printf("!!!!  isnan in candidate point!\n");
            exit(1);
            break;
        }
        if ((candidate[i] < data.lowerbound[i])||(candidate[i] > data.upperbound[i])) return -1;
    }
    return 0;
}

void chaintask(float in_tparam[], int *pdim, int *pnsteps, float *out_tparam, int winfo[4],
        float *init_mean, float *chain_cov)
{
    int i,step;
    int nsteps = *pnsteps;
    int gen_id = winfo[0];
    int chain_id = winfo[1];

    long me = torc_worker_id();

    float leader[data.Nth], loglik_leader;        /* old*/
    float candidate[data.Nth], loglik_candidate;    /* new*/

    for (i = 0; i < data.Nth; i++)
        leader[i] = in_tparam[i];    /*chainwork->in_tparam[i];*/    /* get initial leader */
    loglik_leader = *out_tparam;    /*chainwork->out_tparam[0];*/    /* and its value */

    float pj = runinfo.p[runinfo.Gen];

    int burn_in = data.burn_in;

    for (step = 0; step < nsteps + burn_in; step++) {
        float chain_mean[data.Nth];
        if (step == 0)
            for (i = 0; i < data.Nth; ++i) chain_mean[i] = init_mean[i];
        else
            for (i = 0; i < data.Nth; ++i) chain_mean[i] = leader[i];

#if 0
        int fail = compute_candidate_cov(candidate, chain_mean, chain_cov);
#else
        int fail = compute_candidate(candidate, chain_mean, 1); // I keep this for the moment, for performance reasons
#endif

        if (!fail)
        {
            evaluate_F(candidate, &loglik_candidate, me, gen_id, chain_id, step, 1);    /* this can spawn many tasks*/

            if (data.ifdump && step >= burn_in) torc_update_full_db(candidate, loglik_candidate, NULL, 0, 0);   /* last argument should be 1 if it is a surrogate */

            /* Decide */
            float logprior_candidate = logpriorpdf(candidate, data.Nth);    /* from PanosA */
            float logprior_leader = logpriorpdf(leader, data.Nth);
            float L;
            L = exp((logprior_candidate-logprior_leader)+(loglik_candidate-loglik_leader)*pj);    /* with exp, with log in logpriorpdf and fitfun */

            if (L > 1) L = 1;
            float P = uniformrand(0,1);
            if (P < L) {
                for (i = 0; i < data.Nth; i++) leader[i] = candidate[i];    /* new leader! */
                loglik_leader = loglik_candidate;
                if (step >= burn_in) {
					float logprior_leader = logpriorpdf(leader, data.Nth);
					torc_update_curgen_db(leader, loglik_leader, logprior_leader);
				}
            }
            else {
                /*increase counter or add the leader again in curgen_db*/
                if (step >= burn_in) {
					float logprior_leader = logpriorpdf(leader, data.Nth);
					torc_update_curgen_db(leader, loglik_leader, logprior_leader);
				}

            }
        }
        else {
            /*increase counter or add the leader again in curgen_db*/
            if (step >= burn_in) {
					float logprior_leader = logpriorpdf(leader, data.Nth);
					torc_update_curgen_db(leader, loglik_leader, logprior_leader);
				}
        }
    }

    return;
}

#if 1
int compar_desc(const void* p1, const void* p2)
{
    int dir = +1;   /* -1: ascending order, +1: descending order */
    sort_t *s1 = (sort_t *) p1;
    sort_t *s2 = (sort_t *) p2;

    if (s1->nsel < s2->nsel) return dir;
    if (s1->nsel > s2->nsel) return -dir;
    /*    if (s1->nsel == s2->nsel) return 0;*/
    return 0;
}
#endif

int prepare_newgen(int nchains, cgdbp_t *leaders)
{
    /* process curgen_db -> calculate statitics */
    /* compute probs based on F values */
    /* draw new samples (nchains or user-specified) */
    /* find unique samples: fill the (new) leaders table */
    /* count how many times they appear -> nsteps */
    /* return the new sample size (number of chains) */

    int i, p;
    int newchains; /* = nchains;*/

    int n = curgen_db.entries;

    float *fj = (float *) malloc(n*sizeof(float));
    unsigned int *sel = (unsigned int *) malloc(n*sizeof(sel));

    float **g_x;
    g_x = (float **)malloc(data.Nth*sizeof(float *));
    for (i = 0; i < data.Nth; i++)
        g_x[i] = (float *)malloc(n*sizeof(float));

    {/*start block*/
        float **x = g_x;

        for (p = 0; p < data.Nth; p++) {
            for (i = 0; i < n; i++) {
                x[p][i] = curgen_db.entry[i].point[p];
            }
        }

        float meanx[data.Nth], stdx[data.Nth];
        for (p = 0; p < data.Nth; p++) {
            meanx[p] = compute_mean(x[p], n);
            stdx[p] = compute_std(x[p], n, meanx[p]);
        }

        printf("CURGEN DB (COMPLE) %d\n", runinfo.Gen);
        print_matrix((char *)"means", meanx, data.Nth);
        print_matrix((char *)"std", stdx, data.Nth);
    }/*end block*/

    if (1)
    {/*start block*/
        float **x = g_x;
        int un = 0, unflag, j;

        for (p = 0; p < data.Nth; p++) {
            x[p][un] = curgen_db.entry[0].point[p];    /* un==0*/
        }
        un++;
        for (i = 1; i < n; i++) {
            float xi[data.Nth];
            for (p = 0; p < data.Nth; p++) {
                xi[p] = curgen_db.entry[i].point[p];
            }
            unflag = 1;    /* is this point unique?*/
            for (j = 0; j < un; j++) {    /* check all the previous unique points*/
                int compflag;
                compflag = 1;    /**/
                for (p = 0; p < data.Nth; p++) {
                    if (fabs(xi[p]-x[p][j]) > 1e-6) {
                        /*if (xi[p] != x[p][j]) {*/
                        compflag = 0;    /* they differ*/
                        break;
                    }
                }

                if (compflag == 1) {
                    unflag = 0;    /* not unique, just found it in the unique points table*/
                    break;
                }
            }
            if (unflag) {    /* unique, put it in the table */
                for (p = 0; p < data.Nth; p++) {
                    x[p][un] = xi[p];
                }
                un++;
            }
        }

        runinfo.currentuniques[runinfo.Gen] = un; /*+ 1*/;
        runinfo.acceptance[runinfo.Gen] = (1.0*runinfo.currentuniques[runinfo.Gen])/data.Num[runinfo.Gen]; /* check this*/

        float meanx[data.Nth], stdx[data.Nth];
        for (p = 0; p < data.Nth; p++) {
            meanx[p] = compute_mean(x[p], un);
            stdx[p] = compute_std(x[p], un, meanx[p]);
        }

        printf("CURGEN DB (UNIQUE) %d: [un = %d]\n", runinfo.Gen, un); /* + 1);*/
        print_matrix((char *)"means", meanx, data.Nth);
        print_matrix((char *)"std", stdx, data.Nth);
    } /* end block*/

    {
      float t0 = torc_gettime();
    for (i = 0; i < n; i++) fj[i] = curgen_db.entry[i].F;    /* separate point from F ?*/
    float t1 = torc_gettime();
    calculate_statistics(fj, n, data.Num[runinfo.Gen], runinfo.Gen, sel);
    float t2 = torc_gettime();
    printf("init + calc stats : %lf + %lf = %lf seconds\n", t2-t1, t1-t0, t2-t0);
    }

    newchains = 0;
    for (i = 0; i < n; i++) {
        if (sel[i] != 0) newchains++;
    }

    //sort_t list[n];
    sort_t *list;
    list = calloc(1, n*sizeof(sort_t));
    for (i = 0; i < n; i++) {
        list[i].idx = i;
        list[i].nsel = sel[i];
        list[i].F = curgen_db.entry[i].F;
    }

#if VERBOSE
    printf("Points before\n");
    for (i = 0; i < n; i++) {
        printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
    }
#endif

    qsort(list, n, sizeof(sort_t), compar_desc);

#if VERBOSE
    printf("Points after\n");
    for (i = 0; i < n; i++) {
        printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
    }
#endif

#if 1   /* UPPER THRESHOLD */
      /* peh:check this */
      /* breaking long chains */
    int initial_newchains = newchains;
    int h_threshold = data.MaxChainLength;    /* peh: configuration file + more balanced lengths */
    for (i = 0; i < initial_newchains; i++) {
        if (list[i].nsel > h_threshold) {
            while (list[i].nsel > h_threshold) {
                list[newchains] = list[i];
                list[newchains].nsel = h_threshold;
                list[i].nsel = list[i].nsel - h_threshold;
                newchains++;
            }
        }
    }

    qsort(list, n, sizeof(sort_t), compar_desc);

#if VERBOSE
    printf("Points broken\n");
    for (i = 0; i < n; i++) {
        printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
    }
#endif

#endif

#if 1   /* LOWER THRESHOLD */
    /* single to float step chains */
    int l_threshold = data.MinChainLength;    /* peh: configuration file + more balanced lengths */
    for (i = 0; i < newchains; i++) {
        if ((list[i].nsel > 0)&&(list[i].nsel < l_threshold)) {
            list[i].nsel = l_threshold;
        }
    }

    qsort(list, n, sizeof(sort_t), compar_desc);

#if VERBOSE
    printf("Points advanced\n");
    for (i = 0; i < n; i++) {
        printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
    }
#endif

#endif





    int ldi;    /* leader index*/
    ldi = 0;
    for (i = 0; i < n; i++) {    /* newleader */
        if (list[i].nsel != 0) {
            int idx = list[i].idx;
            for (p = 0; p < data.Nth; p++) {
                leaders[ldi].point[p] = curgen_db.entry[idx].point[p];
            }
            leaders[ldi].F = curgen_db.entry[idx].F;
            leaders[ldi].nsel = list[i].nsel;
            ldi++;
        }
    }

    free(list);

    for (i = 0; i < newchains; i++) leaders[i].queue = -1;    /* rr*/

#if VERBOSE
    printf("Leaders before\n");
    for (i = 0; i < newchains; i++) {
        printf("%d %d %f %d\n", i, leaders[i].nsel, leaders[i].F, leaders[i].queue);
    }
#endif


#if defined(_USE_TORC_)
    /* cool and greedy partitioning ala Panos-- ;-) */

    int nworkers = torc_num_workers();
    int *workload = (int *)calloc(1, nworkers*sizeof(int));    /* workload[1..workers] = 0*/

    for (i = 0; i < newchains; i++) {
        int least_loader_worker = compute_min_idx_i(workload, nworkers);
        leaders[i].queue = least_loader_worker;
        workload[least_loader_worker] += leaders[i].nsel;
    }

    print_matrix_i((char *)"initial workload", workload, nworkers);
    free(workload);
#endif

#if VERBOSE
    printf("Leaders after\n");
    for (i = 0; i < newchains; i++) {
        printf("%d %d %f %d\n", i, leaders[i].nsel, leaders[i].F, leaders[i].queue);
    }
#endif

    if (1)
    {/*start block*/
        /*    float x[data.Nth][n];*/
        float **x = g_x;
        for (i = 0; i < newchains; i++) {
            for (p = 0; p < data.Nth; p++) {
                x[p][i] = leaders[i].point[p];
            }
        }

        float meanx[data.Nth], stdx[data.Nth];
        for (p = 0; p < data.Nth; p++) {
            meanx[p] = compute_mean(x[p], newchains);
            stdx[p] = compute_std(x[p], newchains, meanx[p]);
        }

        printf("CURGEN DB (LEADER) %d: [nlead=%d]\n", runinfo.Gen, newchains);
        print_matrix((char *)"means", meanx, data.Nth);
        print_matrix((char *)"std", stdx, data.Nth);
    }/*end block*/

    if (data.use_local_cov)
        precompute_chain_covariances(leaders, data.init_mean, data.local_cov, newchains);

    curgen_db.entries = 0;    /* reset curgen db*/
    printf("calculate_statistics: newchains=%d\n", newchains);

    for (i = 0; i < data.Nth; i++) free(g_x[i]);
    free(g_x);

    free(fj);
    free(sel);

    return newchains;
}

/*** HELPFUL ***/
void call_gsl_rand_init()
{
    printf("CALLING gsl_rand_init() on node %d\n", torc_node_id()); fflush(0);
    gsl_rand_init(data.seed);
}

void spmd_gsl_rand_init()
{
#if defined(_USE_TORC_)
    int i;
    for (i = 0; i < torc_num_nodes(); i++) {
        torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())call_gsl_rand_init, 0);
    }
    torc_waitall();
#else
    call_gsl_rand_init();
#endif
}

void call_print_matrix_2d()
{
    printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    print_matrix_2d((char *)"runinfo.SS", runinfo.SS, data.Nth, data.Nth);
    printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
}

void spmd_print_matrix_2d()
{
#if defined(_USE_TORC_)
    int i;
    for (i = 0; i < torc_num_nodes(); i++) {
        torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())call_print_matrix_2d, 0);
    }
    torc_waitall();
#else
    call_print_matrix_2d();
#endif
}

void call_update_gdata()    /* step for p[j]*/
{
#if defined(_USE_TORC_)
    MPI_Bcast(runinfo.SS[0], data.Nth*data.Nth, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(runinfo.p, data.MaxStages, MPI_DOUBLE, 0, MPI_COMM_WORLD);    /* just p[Gen]*/
    MPI_Bcast(&runinfo.Gen, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
}

void spmd_update_gdata()    /* step*/
{
#if defined(_USE_TORC_)
    int i;
    if (torc_num_nodes() == 1) return;
    for (i = 0; i < torc_num_nodes(); i++) {
        torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())call_update_gdata, 0);
    }
    torc_waitall();
#endif
}

int main(int argc, char *argv[])
{
    int i;
    float t0, gt0, gt1;
    int winfo[4];
    int nchains; /* was below*/

#if defined(_USE_TORC_)
    torc_register_task((void *)initchaintask);
    torc_register_task((void *)chaintask);
    torc_register_task((void *)torc_update_full_db_task);
    torc_register_task((void *)torc_update_curgen_db_task);
    torc_register_task((void *)torc_update_curres_db_task);
    torc_register_task((void *)reset_nfc_task);
    torc_register_task((void *)get_nfc_task);
    torc_register_task((void *)taskfun);
    torc_register_task((void *)call_gsl_rand_init);
    torc_register_task((void *)call_print_matrix_2d);
    torc_register_task((void *)call_update_gdata);
#endif

    data_init();

    fitfun_initialize(NULL);

#if defined(_USE_TORC_)
    torc_init(argc, argv, MODE_MS);
#endif

#if defined(_AFFINITY_)
    spmd_setaffinity();
#endif

    spmd_gsl_rand_init();

    curgen_db.entries = 0; /* peh+ */

    int goto_next = 0;
    int res;
    if (data.restart) {
    res = load_runinfo();
    if (res == 0) {
        load_curgen_db(runinfo.Gen);
        nchains = data.Num[0];
        printf("nchains = %d\n", nchains);
        gt0 = t0 = torc_gettime();
        goto_next = 1;
    }
    }

    gt0 = t0 = torc_gettime();

    nchains = data.Num[0];
    float out_tparam[data.PopSize];    /* nchains*/

    if (goto_next == 0)
    {

        FILE *init_fp = NULL;
        if (data.load_from_file == 1) {
            init_fp = fopen("init_db.txt", "r");    /* peh: parametrize this */
            if (init_fp == NULL) {
                printf("init_db.txt file not found!\n");
                exit(1);
            }
        }

#if defined(_USE_OPENMP_)
#pragma omp parallel
	{
    printf("Hello from thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());
#pragma omp for
//	{
#endif

        for (i = 0; i < nchains; i++) {
            winfo[0] = runinfo.Gen;
            winfo[1] = i;
            winfo[2] = -1;
            winfo[3] = -1;

            float in_tparam[data.Nth];

            if (data.prior_type == 0)    /* uniform */
            {
                /* uniform */
                int d;
                for (d = 0; d < data.Nth; d++) {
                    in_tparam[d] = uniformrand(0,1);
                    in_tparam[d] *= (data.upperbound[d]-data.lowerbound[d]);
                    in_tparam[d] += data.lowerbound[d];
                }
            }
            else if (data.prior_type == 1)    /* gaussian */
            {
                mvnrnd(data.prior_mu, data.prior_sigma, in_tparam, data.Nth);
            }
            else if (data.prior_type == 3)    /* composite */
            {
                int d;
                for (d = 0; d < data.Nth; d++) {
                    if (data.compositeprior_distr[d] == 0) {
                        in_tparam[d] = uniformrand(0,1);
                        in_tparam[d] *= (data.upperbound[d]-data.lowerbound[d]);
                        in_tparam[d] += data.lowerbound[d];
                    }
                    else if (data.compositeprior_distr[d] == 1) {
                        mvnrnd(&data.prior_mu[d], &data.prior_sigma[d], &in_tparam[d], 1);
                    }
                    else if (data.compositeprior_distr[d] == 2) {
                      	in_tparam[d] = truncated_normal_rand (data.prior_mu[d], data.prior_sigma[d], data.lowerbound[d], data.upperbound[d]);
                    }
                }
            }

            if (data.load_from_file == 1)    /* override the computed point and read it from the file */
            {
              	printf("reading from init_db.txt\n");
                int j;
                for (j = 0; j < data.Nth; j++) fscanf(init_fp, "%lf", &in_tparam[j]);
                fscanf(init_fp, "%lf", &out_tparam[i]);
                float prior_val;
                fscanf(init_fp, "%lf", &prior_val);

                /*torc_update_curgen_db(in_tparam, out_tparam[i]);*/    /* peh - eval or not */
                /*if (data.ifdump) torc_update_full_db(in_tparam, out_tparam[i], NULL, 0, 0);*/
            }

            if (data.prior_type <= 3)    /* peh: file without or with function evaluations? */
            {
#if defined(_USE_TORC_)
                torc_create(-1, (void (*)())initchaintask, 4,
                        data.Nth, MPI_DOUBLE, CALL_BY_COP,
                        1, MPI_INT, CALL_BY_COP,
                        1, MPI_DOUBLE, CALL_BY_RES,
                        4, MPI_INT, CALL_BY_COP,
                        in_tparam, &data.Nth, &out_tparam[i], winfo);
#else

//#if defined(_USE_OPENMP_)
//#endif

//                #pragma omp task shared(data, out_tparam) firstprivate(i, in_tparam, winfo)
//                {
//                initchaintask(
//                        in_tparam, &data.Nth, &out_tparam[i], winfo);
//                }

                #pragma omp task firstprivate(i, winfo, in_tparam) shared(data, out_tparam)
                {
                  //usleep(10*1000);
                  //printf("hello from task %d on worker %d\n", i, omp_get_thread_num());
                  initchaintask(in_tparam, &data.Nth, &out_tparam[i], winfo);
                }

#endif
            }
        }

#if defined(_USE_OPENMP_)
	//} // single
	} // parallel
#endif


#if defined(_USE_TORC_)
        if (data.stealing)
	    torc_enable_stealing();

        torc_waitall();

        if (data.stealing)
            torc_disable_stealing();
#endif

        if (data.load_from_file == 1) {
            fclose(init_fp);
        }

        gt1 = torc_gettime();
        int g_nfeval = get_nfc();
        printf("server: Generation %d: total elapsed time = %lf secs, generation elapsed time = %lf secs for function calls = %d\n", runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
        reset_nfc();

        if (data.icdump) dump_curgen_db(runinfo.Gen);
        if (data.ifdump) dump_full_db(runinfo.Gen);

        /* save here */
        if (data.restart) {
            save_runinfo();
            check_for_exit();
        }

        if (data.MaxStages == 1) goto end;

    }
    ;

    static cgdbp_t *leaders; /*[MAXCHAINS];*/
    leaders = (cgdbp_t *)calloc(1, data.PopSize*sizeof(cgdbp_t));
    for (i = 0; i < data.PopSize; i++) {
        leaders[i].point = (float *)calloc(1, data.Nth*sizeof(float));
    }

    curres_db.entries = 0;
    nchains = prepare_newgen(nchains, leaders);    /* calculate statistics */

    spmd_update_gdata();
    call_print_matrix_2d();

    /* this can be moved above */
    if (runinfo.p[runinfo.Gen] == 1) {
        printf("p == 1 from previous run, nothing more to do\n");
        goto end;
    }

    runinfo.Gen++;
    for (; runinfo.Gen < data.MaxStages; runinfo.Gen++){

        /* process current generation, compute probs, find new chains */
        /*leader[i]: { point[data.Nth], F, nsteps}*/

        int nsteps;
        gt0 = torc_gettime();


#if defined(_USE_OPENMP_)
#pragma omp parallel
	{
#pragma omp single nowait
	{
#endif
        int winfo[4];
        float in_tparam[data.Nth];
        float init_mean[data.Nth];
        float chain_cov[data.Nth*data.Nth];

        for (i = 0; i < nchains; i++) {
            winfo[0] = runinfo.Gen;
            winfo[1] = i;
            winfo[2] = -1;    /* not used */
            winfo[3] = -1;    /* not used */

            int p;
            for (p = 0; p < data.Nth; p++)
                in_tparam[p] = leaders[i].point[p];
            nsteps = leaders[i].nsel;

            if (data.use_local_cov) {
                for (p = 0; p < data.Nth*data.Nth; ++p)
                    chain_cov[p] = data.local_cov[i][p];

                for (p = 0; p < data.Nth; ++p) {
                    if (data.use_proposal_cma)
                        init_mean[p] = data.init_mean[i][p];
                    else
                        init_mean[p] = leaders[i].point[p];
                }
            }
            else {
                int j;
                for (p = 0; p < data.Nth; ++p)
                    for (j = 0; j < data.Nth; ++j)
                        chain_cov[p*data.Nth+j]= data.bbeta*runinfo.SS[p][j];

                for (p = 0; p < data.Nth; ++p)
                    init_mean[p] = in_tparam[p];
            }

            out_tparam[i] = leaders[i].F;    /* loglik_leader...*/

#if defined(_USE_TORC_)
            torc_create(leaders[i].queue, (void (*)())chaintask, 7,
                    data.Nth, MPI_DOUBLE, CALL_BY_COP,
                    1, MPI_INT, CALL_BY_COP,
                    1, MPI_INT, CALL_BY_COP,
                    1, MPI_DOUBLE, CALL_BY_REF,
                    4, MPI_INT, CALL_BY_COP,
                    data.Nth, MPI_DOUBLE, CALL_BY_COP,
                    data.Nth*data.Nth, MPI_DOUBLE, CALL_BY_COP,
                    in_tparam, &data.Nth, &nsteps, &out_tparam[i], winfo,
                    init_mean, chain_cov);
#else

#if defined(_USE_OPENMP_)
	    #pragma omp task shared(data, out_tparam) firstprivate(i, nsteps, in_tparam, winfo, init_mean, chain_cov)
#endif
            chaintask(
                    in_tparam, &data.Nth, &nsteps, &out_tparam[i], winfo,
                    init_mean, chain_cov);
#endif
        }

        /* wait for all chain tasks to finish */

#if defined(_USE_OPENMP_)
	} // single
	} // parallel
#endif


#if defined(_USE_TORC_)
        if (data.stealing)
            torc_enable_stealing();

        torc_waitall();
        if (data.stealing)
            torc_disable_stealing();
#endif

        gt1 = torc_gettime();
        int g_nfeval = get_nfc();
        printf("server: Generation %d: total elapsed time = %lf secs, generation elapsed time = %lf secs for function calls = %d\n", runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
        reset_nfc();

        if (data.icdump) dump_curgen_db(runinfo.Gen);
        if (data.ifdump) dump_full_db(runinfo.Gen);

        /* save here*/
        if (data.restart) {
            save_runinfo();
            check_for_exit();
        }

        curres_db.entries = 0;
        nchains = prepare_newgen(nchains, leaders);    /* calculate statistics*/

        spmd_update_gdata();
        /*spmd_print_matrix_2d();*/
        call_print_matrix_2d();

#if 0
        printf("=================\n");
        print_matrix((char *)"runinfo.p", runinfo.p, runinfo.Gen+1);
        print_matrix((char *)"runinfo.CoefVar", runinfo.CoefVar, runinfo.Gen+1);
        print_matrix_i((char *)"runinfo.currentuniques", runinfo.currentuniques, runinfo.Gen+1);
        print_matrix((char *)"runinfo.acceptance", runinfo.acceptance, runinfo.Gen+1);
        print_matrix((char *)"runinfo.logselection", runinfo.logselection, runinfo.Gen+1);
        printf("=================\n");
#endif

        if (runinfo.p[runinfo.Gen] == 1) {
            break;
        }
        if (runinfo.Gen+1 == data.MaxStages) {
            break;
        }
    }

    if (data.MaxStages == 1) runinfo.Gen = 0;	// small correction for this extreme case

    print_matrix((char *)"runinfo.p", runinfo.p, runinfo.Gen+1);
    print_matrix((char *)"runinfo.CoefVar", runinfo.CoefVar, runinfo.Gen+1);
    print_matrix_i((char *)"runinfo.currentuniques", runinfo.currentuniques, runinfo.Gen+1);
    print_matrix((char *)"runinfo.acceptance", runinfo.acceptance, runinfo.Gen+1);
    print_matrix((char *)"runinfo.logselection", runinfo.logselection, runinfo.Gen+1);

    float logEvidence[1];
    logEvidence[0] = compute_sum(runinfo.logselection, runinfo.Gen+1);
    print_matrix((char *)"logEvidence", logEvidence, 1);

    /* peh:check -- inner tmcmc */
    {
        FILE *fp;
        fp = fopen("fitness.txt", "w");
        fprintf(fp, "%lf\n", logEvidence[0]);
        fclose(fp);
    }

    print_matrix_2d((char *)"runinfo.SS", runinfo.SS, data.Nth, data.Nth);

    for (i = 0; i < runinfo.Gen+1; i++) {
        char title[64];
        sprintf(title, "runinfo.meantheta(%d)", i);
        /*print_matrix((char *)"runinfo.meantheta", runinfo.meantheta[i], data.Nth);*/
        print_matrix((char *)title, runinfo.meantheta[i], data.Nth);
    }

    /* last save here - do we need this? what happens if I restart the program with this saved data*/
    if (data.restart)
        save_runinfo();

end:
    /* making a copy of last curgen_db file */
    if (data.icdump)
    {
        printf("lastgen = %d\n", runinfo.Gen);
        char cmd[256];
        sprintf(cmd, "cp curgen_db_%03d.txt final.txt", runinfo.Gen);
        system(cmd);
    }


    /* shutdown */
    fitfun_finalize();

    printf("total function calls = %d\n", get_tfc());
#if defined(_USE_TORC_)
    torc_finalize();
#endif
    return 0;
}
