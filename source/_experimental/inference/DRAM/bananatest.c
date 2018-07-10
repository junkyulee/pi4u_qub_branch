#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
// testing MCMC with Gaussian target
// high dimensional Gaussian target with positivity constraint

int compare (const void * a, const void * b)
{
	return (*(float*)a - *(float*)b);
}

float sign(float a)
{
	if (a > 0) return +1;
	if (a < 0) return -1;
	return 0;
}

float norm(float *a, int n)
{
	float s = 0;
	for (int i = 0; i < n; i++)
		s += pow(a[i],2.0);

	return sqrt(s);
}

void covcond(float c, float *a, float *Sig, float *Lam, int n)
{
	//%COVCOND covariance matrix with given condition number
	//% [Sig,Lam] = covcond(C,A) generates covariance matrix and its
	//% inverse with given cond number C and first direction A.

	//% create orthogonal basis z, with 1 direction given by 'a'

//	e     = sort(1./linspace(c,1,length(a)));
	float e[n];

	for (int i = 0; i < n; i++)
	{
		e[i] = c + i*((1-c)/(n-1));
		e[i] = 1.0/e[i];
#if DEBUG
		printf("e[%d] = %f\n", i, e[i]);
#endif
	}

	qsort (e, n, sizeof(float), compare);
#if DEBUG
	for (int i = 0; i < n; i++)
	{
		printf("e[%d] = %f\n", i, e[i]);		
	}

	for (int i = 0; i < n; i++)
	{
		printf("a[%d] = %f\n", i, a[i]);		
	}
#endif

//	a(1)  = a(1) + sign(a(1)) * norm(a);  %the Householder trick
	a[0] = a[0] + sign(a[0]) * norm(a, n);

	for (int i = 0; i < n; i++)
	{
		printf("a[%d] = %f\n", i, a[i]);		
	}


//	z     = eye(length(a)) - 2.0/norm(a)^2*a*a';
	float norma2 = pow(norm(a, n), 2);
	float z[n][n];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i == j)
				z[i][j] = 1.0 + (-2.0/norma2)*a[i]*a[j];
			else
				z[i][j] = (-2.0/norma2)*a[i]*a[j];

#if DEBUG
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("z[%d][%d] = %f\n", i, j, z[i][j]);		
	}
#endif

	float S[n][n], L[n][n];
	float ze[n][n], zt[n][n];

//	Sig   = z * diag(e) * z' ;              % target covariance

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ze[i][j] = z[i][j]*e[j];
			zt[i][j] = z[j][i];
		}
	}

#if DEBUG
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("ze[%d][%d] = %f\n", i, j, ze[i][j]);		
	}
#endif

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			float s = 0.0;
			for (int k = 0; k < n; k++)
				s += ze[i][k]*zt[k][j];
			S[i][j] = s;
		}
	}

#if DEBUG
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("S[%d][%d] = %f\n", i, j, S[i][j]);		
	}
#endif

//	Lam   = z * inv(diag(e)) * z';          % and its inverse
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ze[i][j] = z[i][j]*(1.0/e[j]);
			zt[i][j] = z[j][i];
		}
	}

#if DEBUG
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("ze[%d][%d] = %f\n", i, j, ze[i][j]);		
	}
#endif

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			float s = 0.0;
			for (int k = 0; k < n; k++)
				s += ze[i][k]*zt[k][j];
			L[i][j] = s;
		}
	}

#if DEBUG
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("L[%d][%d] = %f\n", i, j, L[i][j]);		
	}
#endif

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Sig[i*n+j] = S[i][j];
			Lam[i*n+j] = L[i][j];
		}
	}
}

//1x2 2x2 = 1x2
//
//[a b]	[1  2]
//	[3  4]  
//[x y] [a b]

float *d_mu;
float *d_imat;
float *d_bpar;

void bananafun(float *y, float *x, float *ab, int inverse)
{
	//%BANANAFUN banana shaped function
	float a = ab[0]; 
	float b = ab[1]; //% parameters

	//y = x;
	y[0] = x[0];
	y[1] = x[1];

	if (inverse)
	{
		y[0] = x[0]/a;
		y[1] = x[1]*a + a*b*(x[0]*x[0] + a*a);
	}
	else
	{
		y[0] = a*x[0];
		y[1] = x[1]/a - b*(y[0]*y[0] + a*a);
	}
}

float ssfun(float *x, int n)
{
	//model.ssfun    = inline('(x-d.mu)*d.Lam*(x-d.mu)''','x','d');

	//bananafun(x-d.mu,d.bpar,1)*d.imat*bananafun(x-d.mu,d.bpar,1)''','x','d');

	float xdmu[n], prod[n];

	for (int i = 0; i < n; i++)
		xdmu[i] = x[i] - d_mu[i];

#if DEBUG
	for (int i = 0; i < n; i++)
		printf("xdmu[%d]=%f\n", i, xdmu[i]);
#endif

	float banana[2];
	bananafun(banana, xdmu, d_bpar, 1);

//1x2 2x2 2x1 
//	( banana[0] banana[1] )  * [	imat[0,0] imat[0,1]
//					imat[1,0] imat[1,1] ]

	for (int i = 0; i < n; i++)
	{
		prod[i] = 0;
		for (int j = 0; j < n; j++)
		{
			prod[i] += banana[j]*d_imat[i*n+j];
		}
	}

#if DEBUG
	for (int i = 0; i < n; i++)
		printf("prod[%d]=%f\n", i, prod[i]);
#endif

	float s = 0;
	for (int i = 0; i < n; i++)
	{
		s += prod[i]*banana[i];
	}
		
#if DEBUG
	printf("s = %f\n", s);
#endif
	return s;
}


float priorfun(float *x, int n)
{
	return 0;
}

struct _params
{
	float *par0;		// = mu+0.1; // initial value
	float sigma2;		// = 1;
	float n0;		// = -1;
	float *lbounds;	// = -inf
	float *ubounds;	// = +inf
	char filename[256];	// = "chain.txt"

	int n;			// = 0;
} params;


void init_params()
{
	params.par0 = NULL;

	params.sigma2 = 1;
	params.n0 = -1;

	params.lbounds = NULL;
	params.ubounds = NULL;
	strcpy(params.filename, "chain.txt");

	params.n = 0;
}

struct _options
{
	int nsimu;		//  =nsimu;
	int adaptint;		// = 500;
	float *qcov;     	//= Sig.*2.4^2./npar;
	float drscale; 	// = drscale;
	float adascale;	// = adascale; // default is 2.4/sqrt(npar) ;

	int printint;		// = 2000;
	float qcovadj;		// = 1e-5;

	int verbosity;		// = 1;
} options;

void init_options()
{
	options.qcovadj = 1e-5;
	options.verbosity = 1;
}


#include "dramrun.c"

int main(int argc, char *argv[])
{
	gsl_rand_init(1);

	//addpath([pwd,filesep,'utils']);

	const int npar     = 2;     // dimension of the target
	float drscale  = 2;     // DR shrink factor
	int adaptint = 10;
	float adascale = 2.4/sqrt(npar); // scale for adaptation
	int nsimu    = 1e5;  // 1e3  // number of simulations
	int pos = 0;           // positivity?


	typedef enum 
	{
		MH,
		DR,
		AM,
		DRAM
	} method_type;

//	method_type method = MH;
//	method_type method = DR;
	method_type method = AM;
//	method_type method = DRAM;

	switch (method)
	{
		case MH: drscale  = 0; adaptint = 0; break;
		case DR: drscale  = 2; adaptint = 0; break;
		case AM: drscale  = 0; adaptint = 10; break;
		case DRAM: drscale  = 2; adaptint = 10; break;
	}

#if 0
	drscale = 2;
	adaptint = -1;
	nsimu = 1e5;
#endif
	printf("drscale = %f, adaptint = %d\n", drscale, adaptint);

	// with positivity, you need nsimu to be at least
	// npar ->  nsimu 
	// 10   ->  10 000
	// 20   ->  500 000
	// 100  ->  5 000 000 ?

	//mu = zeros(1,npar);       // center point
	float *mu = (float *)malloc(npar*sizeof(float));
	for (int i = 0; i < npar; i++)
		mu[i] = 0;

	//cmat = [1 0.9;0.9 1]; % target covariance
	float cmat[npar*npar], imat[npar*npar];
	cmat[0] = 1;
	cmat[1] = 0.9;
	cmat[2] = 0.9;
	cmat[3] = 1;
	inv(imat, cmat, npar);

	//  bpar = [1 1]; % "bananity" of the target, see bananafun.m
	float bpar[2];
	bpar[0] = 1;
	bpar[1] = 1;

	// create input arguments for the dramrun function
	//clear model params options

	// sum of squares function
	//model.ssfun    = inline('(x-d.mu)*d.Lam*(x-d.mu)''','x','d');

	init_params();
	params.par0    = (float *)malloc(npar*sizeof(float));
	for (int i = 0; i < npar; i++)
		params.par0[i] = mu[i];	//+0.1;	// initial value

	params.lbounds = (float *)malloc(npar*sizeof(float));
	for (int i = 0; i < npar; i++)
		params.lbounds[i] = -1e10;	// -inf

	params.ubounds = (float *)malloc(npar*sizeof(float));
	for (int i = 0; i < npar; i++)
		params.ubounds[i] = +1e10;	// -inf

	// positivity
	if (pos)
	{
	//	params.bounds = (ones(npar,2)*diag([0,Inf]));
		for (int i = 0; i < npar; i++)
			params.lbounds[i] = 0;
	}

	strcpy(params.filename, "bananatest_chain.txt");

	// arguments for ssfun are in data
	//data = struct('mu',mu,'imat',imat,'bpar',bpar);
	d_mu = mu;
	d_imat = imat;
	d_bpar = bpar;

	//float x[2] = {0.5, 0.2};
	//float s = ssfun(x, npar);
	//printf("s = %f\n", s);


	init_options();
	options.nsimu    = nsimu;
	options.adaptint = 500;

	//options.qcov     = eye(2)*5;
	options.qcov = (float *)malloc(npar*npar*sizeof(float));
	for (int i = 0; i < npar; i++)
	for (int j = 0; j < npar; j++)
		if (i == j) 
			options.qcov[i*npar+j] = 5.0;
		else
			options.qcov[i*npar+j] = 0.0;

	options.drscale  = drscale;
	options.adaptint = adaptint;
	options.adascale = adascale; // default is 2.4/sqrt(npar) ;
	options.printint = 2000;

	for (int i = 0; i < npar; i++)
		for (int j = 0; j < npar; j++)
			printf("options.qcov[%d][%d]=%f\n", i, j, options.qcov[i*npar+j]);



//	{
//	float xpar[2] = {0.5, 0.5};
//	float s = ssfun(xpar, npar);
//	printf("ssfun(%f,%f) =  %f\n", xpar[0], xpar[1], s);
//	}

	dramrun(npar);

//	[results,chain] = dramrun(model,data,params,options);


	return 0;
}
