/*
 *  engine_stats.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <math.h>
#include "engine_tmcmc.h"
#include <time.h>

int display = 0;

/** OBJLOGP FUNCTION **/
float Objlogp(float x, float *fj, int fn, float pj, float tol)
{
	int i;
	float fjmax = compute_max(fj, fn);

	float weight[fn];
	for (i = 0; i < fn; i++)
		weight[i] = exp((fj[i]-fjmax)*(x-pj));

	float sum_weight = compute_sum(weight, fn);

	float q[fn];
	for (i = 0; i < fn; i++)
		q[i] = weight[i]/sum_weight;

	float mean_q = compute_mean(q, fn);
	float std_q = compute_std(q, fn, mean_q);

	float CoefVar = pow(std_q/mean_q-tol, 2);	// result

	return CoefVar;
}

typedef struct fparam_s {
	float *fj;
	int     fn;
	float  pj;
	float  tol;
} fparam_t;

//fparam_t *sfp;

//float Objlogp_s(float *x, int n)
//{
//	float *fj = sfp->fj;
//	int fn = sfp->fn;
//	float pj = sfp->pj;
//	float tol = sfp->tol;
//
//	return Objlogp(x[0], fj, fn, pj, tol);
//}

float Objlogp_gsl(float x, void *param)
{
	fparam_t *fp = (fparam_t *) param;

	float *fj = fp->fj;
	int fn = fp->fn;
	float pj = fp->pj;
	float tol = fp->tol;

	float res = Objlogp(x, fj, fn, pj, tol);
//	printf("Objlogp(%lf)=%lf\n", x, res);
	return res;
}

float Objlogp_gsl2(const gsl_vector *v, void *param)
{
	float x;
	x = gsl_vector_get(v, 0);

	return Objlogp_gsl(x, param);
}

/*** OPTIMIZATION ***/
int fminsearch(tmcmc_data_t *tmcmc_data, float *fj, int fn, float pj, float tol, float *xmin, float *fmin)
{
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	int conv = 0;

	size_t iter = 0, max_iter = tmcmc_data->data.options.MaxIter;	// USER input
	float Tol = tmcmc_data->data.options.Tol;
	int Display = tmcmc_data->data.options.Display;
	int status;
	float size;

	fparam_t fp;

	fp.fj = fj; fp.fn = fn; fp.pj = pj; fp.tol = tol;

	/* Starting point */
	x = gsl_vector_alloc (1);
	gsl_vector_set (x, 0, pj);

	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (1);
	gsl_vector_set_all (ss, 0.001);	// input?

       /* Initialize method and iterate */
	minex_func.n = 1;
	minex_func.f = Objlogp_gsl2;
	minex_func.params = &fp;

       s = gsl_multimin_fminimizer_alloc (T, 1);
       gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, Tol);

		if (status == GSL_SUCCESS) {
			conv = 1;
			if (Display)
				printf ("converged to minimum at\n");
		}

		if (Display)
			printf ("%3ld %.16lf f() = %.16f size = %.16f\n",
				iter, gsl_vector_get (s->x, 0), s->fval, size);

	} while (status == GSL_CONTINUE && iter < max_iter);

	if (conv) {
		*fmin = s->fval;
		*xmin = gsl_vector_get(s->x, 0);
	} else {
		*fmin = 0;
		*xmin = 0.0;
	}

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return conv;
}

int fmincon(tmcmc_data_t *tmcmc_data, float *fj, int fn, float pj, float tol, float *xmin, float *fmin)
{
	int status;
	int iter = 0, max_iter = tmcmc_data->data.options.MaxIter;	// USER input
	float Tol = tmcmc_data->data.options.Tol;
	int Display = tmcmc_data->data.options.Display;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	float m = 0.5;     // input
 	float a = 0.0, b = 2.0;    // input
	gsl_function F;
	int conv = 0;
	gsl_vector *x;
	int i;

	x = gsl_vector_alloc (1);

	fparam_t fp;

	fp.fj = fj; fp.fn = fn; fp.pj = pj; fp.tol = tol;

	F.function = Objlogp_gsl;
	F.params = &fp;

	T = gsl_min_fminimizer_goldensection; //brent;
	s = gsl_min_fminimizer_alloc (T);

//	printf("f(a)=%lf\n", Objlogp_gsl(a, &fp));
//	printf("f(b)=%lf\n", Objlogp_gsl(b, &fp));
//	printf("f(m)=%lf\n", Objlogp_gsl(m, &fp));
	float fa = Objlogp_gsl(a, &fp);
	float fb = Objlogp_gsl(b, &fp);
	for (i = 0; i < max_iter; i++) {
		m = a + i*(b-a)/max_iter;
		float fm = Objlogp_gsl(m, &fp);
		if ((fm < fa) && (fm < fb)) break;
	}

	if (i == max_iter) {
		if (Display)
			printf("failed to initialize fmincon!\n");
		return 0;
	}
	else {
		if (Display)
			printf("inited with %d tries\n", i);
	}

	gsl_min_fminimizer_set (s, &F, m, a, b);
//	printf("bbb\n");

	if (Display) {
		printf ("using %s method\n", gsl_min_fminimizer_name (s));
		printf ("%5s [%18s, %18s] %18s %18s\n", "iter", "lower", "upper", "min", "err(est)");
		printf ("%5d [%.16f, %.16f] %.16f %.16f\n", iter, a, b, m, b - a);
	}

	do {
		iter++;
		status = gsl_min_fminimizer_iterate (s);

		m = gsl_min_fminimizer_x_minimum (s);
		a = gsl_min_fminimizer_x_lower (s);
		b = gsl_min_fminimizer_x_upper (s);

		status = gsl_min_test_interval (a, b, Tol, 0.0);
		if (status == GSL_SUCCESS) {
			if (Display)
				printf ("Converged:\n");
			conv = 1;
		}

		if (Display)
			printf ("%5d [%.16f, %.16f]  %.16f %.16f\n",
				iter, a, b, m, b - a);

	} while (status == GSL_CONTINUE && iter < max_iter);

	if (conv) {
		gsl_vector_set (x, 0, m);
		*fmin = Objlogp_gsl(m, &fp);
		*xmin = m;
	} else {
		*fmin = 0;
		*xmin = 0.0;
	}

	gsl_vector_free(x);
	gsl_min_fminimizer_free (s);

	return conv;

}


/*** STATISTICS ***/
void calculate_statistics(tmcmc_data_t *tmcmc_data, float flc[], int n, int nselections, int gen, unsigned int sel[])
{
	//float pflag = 0;
	float tolCOV = tmcmc_data->data.TolCOV;
	float *CoefVar = tmcmc_data->runinfo.CoefVar;
	float *p = tmcmc_data->runinfo.p;
	int *Num = tmcmc_data->data.Num;
	int *currentuniques = tmcmc_data->runinfo.currentuniques;
	float *logselection = tmcmc_data->runinfo.logselection;

	float fmin = 0, xmin = 0;
	int conv = 0;
//	conv = fmincon(flc, n, p[gen], tolCOV, &xmin, &fmin);
	if (display)
		printf("fmincon: conv=%d xmin=%.16lf fmin=%.16lf\n", conv, xmin, fmin);
	if (!conv) {
		conv = fminsearch(tmcmc_data, flc, n, p[gen], tolCOV, &xmin, &fmin);
		if (display)
			printf("fminsearch: conv=%d xmin=%.16lf fmin=%.16lf\n", conv, xmin, fmin);
	}

	/* gen: next generation number */
	int j = gen+1;

	p[j] = xmin;
	CoefVar[j] = fmin;

	if (p[j]>1) {
		//pflag=p[j-1];
		p[j] = 1;
		Num[j]=tmcmc_data->data.LastNum;
	}

	//print_matrix("p", p, j);

	//if (p[j]>10) { // tmcmc_data->data.pstrict
		// tmcmc_data->data.Rsq = tmcmc_data->data.Rsqstrict;
	//}

	// Compute weights and normalize
	float fjmax = compute_max(flc, n);

	int i;

	float weight[n];
	for (i = 0; i < n; i++)
		weight[i] = exp((flc[i]-fjmax)*(p[j]-p[j-1]));

	if (display)
	print_matrix("weight", weight, n);

	float sum_weight = compute_sum(weight, n);

	float q[n];
	for (i = 0; i < n; i++)
		q[i] = weight[i]/sum_weight;

	if (display)
		print_matrix("runinfo_q", q, n);

	//float sum_q = compute_sum(q, n);

	logselection[gen] = log(sum_weight/currentuniques[gen])+fjmax*(p[gen+1]-p[gen]);

	if (display)
		print_matrix("logselection", logselection, gen+1);

	float mean_q = compute_mean(q, n);
	float std_q = compute_std(q, n, mean_q);

	CoefVar[gen] = std_q/mean_q;

	if (display)
		print_matrix("CoefVar", CoefVar, gen+1);

	size_t K = n;
	unsigned int N = 1;

	unsigned int samples = n; //1000;
	unsigned int nn[samples];

	for (i = 0; i < (int)samples; i++) sel[i] = 0;

	int k;

	if (nselections == 0) nselections = samples; // n;
	for (k = 0; k < nselections; k++) {

		//gsl_ran_multinomial (r, K, N, q, nn);
		multinomialrand (K, N, q, nn);
		for (i = 0; i < (int)K; i++) sel[i]+=nn[i];
	}

	if (display) {
		printf("\n s = [");
		for (i = 0; i < (int)K; i++) printf("%d ", sel[i]);
		printf("]\n");
	}

	/* compute SS */
	int PROBDIM = tmcmc_data->data.Nth;

	float mean_of_theta[PROBDIM];

	for (i = 0; i < PROBDIM; i++) {
		mean_of_theta[i] = 0;
		for (j = 0; j < n; j++) mean_of_theta[i]+=tmcmc_data->curgen_db.entry[j].point[i]*q[j];

		tmcmc_data->runinfo.meantheta[gen][i] = mean_of_theta[i];
	}

	if (display)
		print_matrix("mean_of_theta", mean_of_theta, PROBDIM);

	float meanv[PROBDIM];
	for (i = 0; i < PROBDIM; i++) {
		meanv[i] = mean_of_theta[i];
	}

	for (i = 0; i < PROBDIM; i++) {
		for (j = 0; j < PROBDIM; j++) {
			float s;
			int k;
			s = 0;
			for (k = 0; k < n; k++) {
				s += q[k]*(tmcmc_data->curgen_db.entry[k].point[i]-meanv[i])*(tmcmc_data->curgen_db.entry[k].point[j]-meanv[j]);
			}
			tmcmc_data->runinfo.SS[i][j] = tmcmc_data->runinfo.SS[j][i] = s;
		}
	}

	if (display)
		print_matrix_2d("runinfo.SS", tmcmc_data->runinfo.SS, PROBDIM, PROBDIM);
}

float priorpdf(float *theta, int n, float *lowerbound, float *upperbound)
{
	int i;
	float res = 1;

	for (i = 0; i < n; i++) {
//		res *= gsl_ran_flat_pdf(theta[i], tmcmc_data.data.lowerbound[i], tmcmc_data.data.upperbound[i]);
		res *= gsl_ran_flat_pdf(theta[i], lowerbound[i], upperbound[i]);
	}
	return res;
}

float posterior(float *theta, int n, float LH)
{
	float res;
//	float Prior = priorpdf(theta, n);
//	printf("Prior = %lf\n", Prior);
//
//	if (Prior > 0)
//		res = LH + log(Prior);	/* xxx */
//	else
//		res = LH;

	res = LH;	/* Algorithm fix by PanosA */

	return res;
}
