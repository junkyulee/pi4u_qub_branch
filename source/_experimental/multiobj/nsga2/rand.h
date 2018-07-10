/* Declaration for random number related variables and routines */

# ifndef _RAND_H_
# define _RAND_H_

/* Variable declarations for the random number generator */
extern float seed;
extern float oldrand[55];
extern int jrand;

/* Function declarations for the random number generator */
void randomize(void);
void warmup_random (float seed);
void advance_random (void);
float randomperc(void);
int rnd (int low, int high);
float rndreal (float low, float high);

# endif
