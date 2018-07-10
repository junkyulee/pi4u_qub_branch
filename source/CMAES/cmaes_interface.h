/* --------------------------------------------------------- */
/* --- File: cmaes_interface.h - Author: Nikolaus Hansen --- */
/* ---------------------- last modified:  IV 2007        --- */
/* --------------------------------- by: Nikolaus Hansen --- */
/* --------------------------------------------------------- */
/*   
     CMA-ES for non-linear function minimization. 

     Copyright (C) 1996, 2003, 2007 Nikolaus Hansen. 
     e-mail: hansen AT lri.fr
     
     Documentation: see file docfunctions.txt
     
     License: see file cmaes.c
*/
#include "cmaes.h"

/* --------------------------------------------------------- */
/* ------------------ Interface ---------------------------- */
/* --------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

/* --- initialization, constructors, destructors --- */
float * cmaes_init(cmaes_t *, int dimension , float *xstart, 
        float *stddev, long seed, int lambda, 
        const char *input_parameter_filename);
void cmaes_init_para(cmaes_t *, int dimension , float *xstart, 
        float *stddev, long seed, int lambda, 
        const char *input_parameter_filename);
float * cmaes_init_final(cmaes_t *);
void cmaes_resume_distribution(cmaes_t *evo_ptr, char *filename);
void cmaes_exit(cmaes_t *);

/* --- core functions --- */
float * const * cmaes_SamplePopulation(cmaes_t *);
float *         cmaes_UpdateDistribution(cmaes_t *, 
                 const float *rgFitnessValues);
const char *     cmaes_TestForTermination(cmaes_t *);

/* --- additional functions --- */
float * const * cmaes_ReSampleSingle( cmaes_t *t, int index);
float const *   cmaes_ReSampleSingle_old(cmaes_t *, float *rgx); 
float *         cmaes_SampleSingleInto( cmaes_t *t, float *rgx);
void             cmaes_UpdateEigensystem(cmaes_t *, int flgforce);

/* --- getter functions --- */
float         cmaes_Get(cmaes_t *, char const *keyword);
const float * cmaes_GetPtr(cmaes_t *, char const *keyword); /* e.g. "xbestever" */
float *       cmaes_GetNew( cmaes_t *t, char const *keyword); /* user is responsible to free */
float *       cmaes_GetInto( cmaes_t *t, char const *keyword, float *mem); /* allocs if mem==NULL, user is responsible to free */

/* --- online control and output --- */
void           cmaes_ReadSignals(cmaes_t *, char const *filename);
void           cmaes_WriteToFile(cmaes_t *, const char *szKeyWord,
                                 const char *output_filename); 
char *         cmaes_SayHello(cmaes_t *);
/* --- misc --- */
float *       cmaes_NewDouble(int n); /* user is responsible to free */
void           cmaes_FATAL(char const *s1, char const *s2, char const *s3, 
               char const *s4);

#ifdef __cplusplus
} // end extern "C"
#endif
