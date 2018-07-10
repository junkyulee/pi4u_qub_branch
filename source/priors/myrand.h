#ifndef _MYRAND_H_
#define _MYRAND_H_

void gsl_rand_init(int seed);

float normal_pdf(float x, float *p);
float normal_log_pdf(float x, float *p);
float normal_rnd( float *p );


float uniform_pdf(float x, float *p);
float uniform_log_pdf(float x, float *p);
float uniform_rnd( float *p );

float exp_pdf(float x, float *p);
float exp_log_pdf(float x, float *p);
float exp_rnd( float *p );

float gamma_pdf(float x, float *p);
float gamma_log_pdf(float x, float *p);
float gamma_rnd( float *p );





float normalrand(float mu, float sigma);
float uniformrand(float a, float b);



#endif

