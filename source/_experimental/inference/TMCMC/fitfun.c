#include <math.h>
#include <unistd.h>
// activate one of the following options
#define _USE_ROSENBROCK_
//#define _USE_BVNPDF_
//#define _USE_MIXED_BVNPDF_
//#define _USE_MIXED_MVNPDF_

#if defined(_USE_MIXED_BVNPDF_)
#include "gsl_headers.h"
float mixedbvnpdf(float *x, int n) /* bivariate */
{
	float P;

	P = gsl_ran_bivariate_gaussian_pdf(x[0]-(-5), x[1]-(-5), 1, 1, 0);
	P += gsl_ran_bivariate_gaussian_pdf(x[0]-(+5), x[1]-(+5), 1, 1, 0);
	return P;
}
#endif

#if defined(_USE_MIXED_MVNPDF_)
#include "gsl_headers.h"
extern float mvnpdf(int n, float *xv, float *mv, float *vm);
float mixedmvnpdf(float *x, int n) /* multivariate */
{
	float P = 0;
	float m1[n];
	float m2[n];

	int i;
	for (i = 0; i < n; i++) m1[i] = -5;
	for (i = 0; i < n; i++) m2[i] = +5;

	P =  mvnpdf(n, x, m1, NULL);
	P += mvnpdf(n, x, m2, NULL);
	return P;
}
#endif


#if defined(_USE_BVNPDF_)
#include "gsl_headers.h"
float bvnpdf(float *x, int n) /* bivariate */
{
	float P;

	P = gsl_ran_bivariate_gaussian_pdf(x[0], x[1], 1, 1, 0);
	return P;
}
#endif

//#include "engine_tmcmc.h"
//extern data_t data;

void fitfun_initialize(char *name)
{
}

void fitfun_finalize()
{
}

float fitfun(float /*const*/ *x, int N, void *output, int *info)
{
	float f;

#if defined(_USE_ROSENBROCK_)
	int i;
	f = 0.0;
	for (i=0; i<N-1; i++)	/* rosenbrock */
		f = f + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
	f = -f;
//	f = -log(f);
#endif

#if defined(_USE_BVNPDF_)
	f = log(bvnpdf(x, N));
#endif

#if defined(_USE_MIXED_BVNPDF_)
	f = log(mixedbvnpdf(x, N));
#endif

#if defined(_USE_MIXED_MVNPDF_)
	f = log(mixedmvnpdf(x, N));
#endif

#if DEBUG
	usleep(100*1000);
#endif
	return f;
}
