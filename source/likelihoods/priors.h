#ifndef _PRIORS_H_
#define _PRIORS_H_

static inline float unifpdf(float x, float l, float u);
static inline float unifpdf2(float x, float l, float len);
static inline float trnpdf(float x, float m, float s, float l, float u);
static inline float log_lognpdf(float x, float m, float s);

//inline float log_normal_pdf(float x, float mu, float sigma);


float priorHB(float *x, float *psi, int n);
float log_priorHB(float *x, float *psi, int n);

#endif
