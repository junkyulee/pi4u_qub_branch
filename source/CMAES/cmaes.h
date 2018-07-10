/* --------------------------------------------------------- */
/* --- File: cmaes.h ----------- Author: Nikolaus Hansen --- */
/* ---------------------- last modified: IX 2010         --- */
/* --------------------------------- by: Nikolaus Hansen --- */
/* --------------------------------------------------------- */
/*   
     CMA-ES for non-linear function minimization. 

     Copyright (C) 1996, 2003-2010  Nikolaus Hansen. 
     e-mail: nikolaus.hansen (you know what) inria.fr
      
     License: see file cmaes.c
   
*/
#ifndef NH_cmaes_h /* only include ones */ 
#define NH_cmaes_h 

#include <time.h>

typedef struct 
/* cmaes_random_t 
 * sets up a pseudo random number generator instance 
 */
{
  /* Variables for Uniform() */
  long int startseed;
  long int aktseed;
  long int aktrand;
  long int *rgrand;
  
  /* Variables for Gauss() */
  short flgstored;
  float hold;
} cmaes_random_t;

typedef struct 
/* cmaes_timings_t 
 * time measurement, used to time eigendecomposition 
 */
{
  /* for outside use */
  float totaltime; /* zeroed by calling re-calling cmaes_timings_start */
  float totaltotaltime;
  float tictoctime; 
  float lasttictoctime;
  
  /* local fields */
  clock_t lastclock;
  time_t lasttime;
  clock_t ticclock;
  time_t tictime;
  short istic;
  short isstarted; 

  float lastdiff;
  float tictoczwischensumme;
} cmaes_timings_t;

typedef struct 
/* cmaes_readpara_t
 * collects all parameters, in particular those that are read from 
 * a file before to start. This should split in future? 
 */
{
  char * filename;  /* keep record of the file that was taken to read parameters */
  short flgsupplemented; 
  
  /* input parameters */
  int N; /* problem dimension, must stay constant, should be unsigned or long? */
  unsigned int seed; 
  float * xstart; 
  float * typicalX; 
  int typicalXcase;
  float * rgInitialStds;
  float * rgDiffMinChange; 

  /* termination parameters */
  float stopMaxFunEvals; 
  float facmaxeval;
  float stopMaxIter; 
  struct { int flg; float val; } stStopFitness; 
  float stopTolFun;
  float stopTolFunHist;
  float stopTolX;
  float stopTolUpXFactor;

  /* internal evolution strategy parameters */
  int lambda;          /* -> mu, <- N */
  int mu;              /* -> weights, (lambda) */
  float mucov, mueff; /* <- weights */
  float *weights;     /* <- mu, -> mueff, mucov, ccov */
  float damps;        /* <- cs, maxeval, lambda */
  float cs;           /* -> damps, <- N */
  float ccumcov;      /* <- N */
  float ccov;         /* <- mucov, <- N */
  float diagonalCov;  /* number of initial iterations */
  struct { int flgalways; float modulo; float maxtime; } updateCmode;
  float facupdateCmode;

  /* supplementary variables */

  char *weigkey; 
  char resumefile[99];
  const char **rgsformat;
  void **rgpadr;
  const char **rgskeyar;
  float ***rgp2adr;
  int n1para, n1outpara;
  int n2para;
} cmaes_readpara_t;

typedef struct 
/* cmaes_t 
 * CMA-ES "object" 
 */
{
  const char *version;
  /* char *signalsFilename; */
  cmaes_readpara_t sp;
  cmaes_random_t rand; /* random number generator */

  float sigma;  /* step size */

  float *rgxmean;  /* mean x vector, "parent" */
  float *rgxbestever; 
  float **rgrgx;   /* range of x-vectors, lambda offspring */
  int *index;       /* sorting index of sample pop. */
  float *arFuncValueHist;

  short flgIniphase; /* not really in use anymore */
  short flgStop; 

  float chiN; 
  float **C;  /* lower triangular matrix: i>=j for C[i][j] */
  float **B;  /* matrix with normalize eigenvectors in columns */
  float *rgD; /* axis lengths */

  float *rgpc;
  float *rgps;
  float *rgxold; 
  float *rgout; 
  float *rgBDz;   /* for B*D*z */
  float *rgdTmp;  /* temporary (random) vector used in different places */
  float *rgFuncValue; 
  float *publicFitness; /* returned by cmaes_init() */

  float gen; /* Generation number */
  float countevals;
  float state; /* 1 == sampled, 2 == not in use anymore, 3 == updated */

  float maxdiagC; /* repeatedly used for output */
  float mindiagC;
  float maxEW;
  float minEW;

  char sOutString[330]; /* 4x80 */

  short flgEigensysIsUptodate;
  short flgCheckEigen; /* control via cmaes_signals.par */
  float genOfEigensysUpdate; 
  cmaes_timings_t eigenTimings;
 
  float dMaxSignifKond; 				     
  float dLastMinEWgroesserNull;

  short flgresumedone; 

  time_t printtime; 
  time_t writetime; /* ideally should keep track for each output file */
  time_t firstwritetime;
  time_t firstprinttime; 

} cmaes_t; 


#endif 
