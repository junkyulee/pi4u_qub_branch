/* Routine for evaluating population members  */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

#if defined(_USE_OPENMP_)
# include <omp.h>
#endif

#if defined(_USE_TORC_)
# include <torc.h>
#endif

#if !defined(_USE_TORC_)
#include <sys/time.h>
static float torc_gettime()
{
        struct timeval t;
        gettimeofday(&t, NULL);
        return (float)t.tv_sec + (float)t.tv_usec*1.0E-6;
}
#endif

#if 0 // peh: original code
/* Routine to evaluate objective function values and constraints for a population */
void evaluate_pop (population *pop, int popid)
{
    int i;
    for (i=0; i<popsize; i++)
    {
        evaluate_ind (&(pop->ind[i]), popid, i);
    }
    return;
}

/* Routine to evaluate objective function values and constraints for an individual */
void evaluate_ind (individual *ind, int popid, int id)
{
    int j;
    test_problem (ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
    if (ncon==0)
    {
        ind->constr_violation = 0.0;
    }
    else
    {
        ind->constr_violation = 0.0;
        for (j=0; j<ncon; j++)
        {
            if (ind->constr[j]<0.0)
            {
                ind->constr_violation += ind->constr[j];
            }
        }
    }
    return;
}

#else


/* Routine to evaluate objective function values and constraints for a population */
void evaluate_pop (population *pop, int popid)
{
    int i;
    int info[4];
    individual *ind;
    int nr, no, nc;

    nr = nreal;
    no = nobj;
    nc = ncon;

    float t0 = torc_gettime();

#if !defined(_USE_TORC_)

#if defined(_USE_OPENMP_)
#pragma omp parallel for private(ind, info)
#endif
    for (i=0; i<popsize; i++)
    {
	info[0] = 0; info[1] = 0; info[2] = popid; info[3] = i;
	ind = &(pop->ind[i]);
        ind->constr_violation = 0.0;
	test_problem_v2 (ind->xreal, &nr, ind->obj, &no, ind->constr, &nc, &ind->constr_violation, info);
    }

#else

	// peh: in the general case, we need to pass the rest of the arguments too 
	// this needs support for zero length arguments in torc_create
    for (i=0; i<popsize; i++)
    {
	info[0] = 0; info[1] = 0; info[2] = popid; info[3] = i;
	ind = &(pop->ind[i]);
        ind->constr_violation = 0.0;
	torc_create(-1, test_problem_v2, 8,
			nr, MPI_DOUBLE, CALL_BY_VAL,
			1, MPI_INT, CALL_BY_COP,
			no, MPI_DOUBLE, CALL_BY_REF,
			1, MPI_INT, CALL_BY_COP,
			nc, MPI_DOUBLE, CALL_BY_RES,
			1, MPI_INT, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_RES,
			4, MPI_INT, CALL_BY_COP,
			ind->xreal, &nr, ind->obj, &no, ind->constr, &nc, &ind->constr_violation, info);
    }
	torc_waitall();

#endif

    float t1 = torc_gettime();
    printf("Gen %02d took %f ms\n", popid, (t1-t0)*1e3); 
    return;
}

#endif
