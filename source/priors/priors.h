#ifndef _PRIORS_H_
#define _PRIORS_H_



typedef struct{
	float (*f)(float, float *);		// Density function
	float (*lf)(float, float *);		// Log of density
	float (*r)( float *);				// Random number
	float *par;						// Parameters of density
	int npar;							// Number of parameters
	char   name[24];					// Name of the distribution
} Density;


//TODO: define the struct Prior:
//typedef struct{
//	Density *d;
//	int Nd;
//} Prior;


void delete_density( Density *d );
void delete_prior( Density *d, int N);

float eval_density(Density d, float x);
float eval_log_density(Density d, float x);
float eval_random( Density d );

void  print_density( Density d );
void  print_priors( Density *d, int N);

float prior_pdf( Density *d, int N, float *x);
float prior_log_pdf( Density *d, int N, float *x);

void read_priors(const char *file,  Density **p_priors, int *p_N );

void reassign_prior( Density *p, int Np, float *psi );
void new_prior_from_prior( Density **new_prior, Density *from_prior, int Npr );

void check_n( int N );



#endif
