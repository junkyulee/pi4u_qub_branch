/*
 *  sa_deriv.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include "fitfun.c" 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <torc.h>
#include "posdef.c"
#include <unistd.h>
int sqrtm(float **A, float **sA, int n);
float* allocate1D(int n);
float** allocate2D(int n);
void fprint_matrix_1d(FILE *fp, char *title, float *v, int n);
void fprint_matrix_2d(FILE *fp, char *title, float **v, int n1, int n2);
void reset_nfc();
void inc_nfc();
int get_nfc();
int get_tfc();
void aux_init();

struct data_s {
	/* read dimension, xl, xu, theta, iord, h0, stepratio, maxsteps  */
	int Nth;
	int iord;

//	float h0;
//	float stepratio;
//	int maxsteps;

	float lb, ub;
	float *lowerbound;
	float *upperbound;
	float *theta;

	float uh;
	int seed;
	int reps;	// total random vectors generated

	float *diffstep;
	int posdef;

} data;

void data_init()
{
	int i;

	data.Nth = 4;

	data.iord = 2;

	data.lb = -1e38;
	data.ub = +1e38;

	data.lowerbound = malloc(data.Nth*sizeof(float));
	data.upperbound = malloc(data.Nth*sizeof(float));
	data.theta = malloc(data.Nth*sizeof(float));

	for (i = 0; i < data.Nth; i++) data.lowerbound[i] = data.lb;
	for (i = 0; i < data.Nth; i++) data.upperbound[i] = data.ub;
	for (i = 0; i < data.Nth; i++) data.theta[i] = 1.0;

	data.uh = 1e-4;
	data.diffstep = malloc(data.Nth*sizeof(float));
	for (i = 0; i < data.Nth; i++) data.diffstep[i] = data.uh;
	data.posdef = -1;
	data.seed = 280675;

//	data.gH_avg = 5000;	

/**********************************************/
// PARAMS
/**********************************************/

	/* USER-DEFINED VALUES */
	FILE *f = fopen("grad.par", "r");
	if (f == NULL) {
		return;
	}

	/*
	Nth		4
	order		2
	Bdef		-4	4

	lb		0 0 0 0 
	ub		2 2 2 2
	theta		1 1 1 1

	Hdef		1e-4
	diffstep	1-e4 1e-4 1e-4 1e-4 
	posdef		-1
	*/

	char line[256];

	int line_no = 0;
	while (fgets(line, 256, f)!= NULL) {
		line_no++;
		if ((line[0] == '#')||(strlen(line)==0)) {
			printf("ignoring line %d\n", line_no);
			continue;
		}

		if (strstr(line, "Nth")) {
			sscanf(line, "%*s %d", &data.Nth);
		}
		else if (strstr(line, "order")) {
			sscanf(line, "%*s %d", &data.iord);
		}
		else if (strstr(line, "Bdef")) {
			sscanf(line, "%*s %lf %lf", &data.lb, &data.ub);
		}
		else if (strstr(line, "Hdef")) {
			sscanf(line, "%*s %lf", &data.uh);
		}
		else if (strstr(line, "posdef")) {
			sscanf(line, "%*s %d", &data.posdef);
		}
		else if (strstr(line, "reps")) {
			sscanf(line, "%*s %d", &data.reps);
		}
		else if (strstr(line, "seed")) {
			sscanf(line, "%*s %d", &data.seed);
		}
	}

	rewind(f);
	line_no = 0;

	free(data.lowerbound);
	free(data.upperbound);
	free(data.theta);
	free(data.diffstep);
	data.lowerbound = (float *)malloc(data.Nth*sizeof(float));
	data.upperbound = (float *)malloc(data.Nth*sizeof(float));
	data.theta = (float *)malloc(data.Nth*sizeof(float));
	data.diffstep = (float *)malloc(data.Nth*sizeof(float));

	for (i = 0; i < data.Nth; i++) data.lowerbound[i] = data.lb;
	for (i = 0; i < data.Nth; i++) data.upperbound[i] = data.ub;
	for (i = 0; i < data.Nth; i++) data.theta[i] = 1.0;
	for (i = 0; i < data.Nth; i++) data.diffstep[i] = data.uh;

	rewind(f);
	while (fgets(line, 256, f)!= NULL) {
		line_no++;

		if ((line[0] == '#')||(strlen(line)==0)) continue;

		if (strstr(line, "theta")) {
			char *p;
			p = strtok(line, " \t\n");
			for (i = 0; i < data.Nth; i++) {
				p = strtok(NULL, " \t\n");
				if (!p) break;
				data.theta[i] = atof(p);
                        }
                }
		else if (strstr(line, "lb")) {
			char *p;
			p = strtok(line, " \t\n");
			for (i = 0; i < data.Nth; i++) {
				p = strtok(NULL, " \t\n");
				if (!p) break;
				data.lowerbound[i] = atof(p);
                        }
                }
		else if (strstr(line, "ub")) {
			char *p;
			p = strtok(line, " \t\n");
			for (i = 0; i < data.Nth; i++) {
				p = strtok(NULL, " \t\n");
				if (!p) break;
				data.upperbound[i] = atof(p);
                        }
                }
		else if (strstr(line, "diffstep")) {
			char *p;
			p = strtok(line, " \t\n");
			for (i = 0; i < data.Nth; i++) {
				p = strtok(NULL, " \t\n");
				if (!p) break;
				data.diffstep[i] = atof(p);
                        }
                }
        }

	for (i = 0; i < data.Nth; i++) {
                printf("param %d: %f < %f < %f: %f\n", i, data.lowerbound[i], data.theta[i], data.upperbound[i], data.diffstep[i]);
        }

	fclose(f);


}

#define NOISE	0

//Noise=@(x) 0.0*randn(size(x));
void Noise(float *x, float *n, int p, unsigned short *dr48_buffer)
{
	int i;
	for (i = 0; i < p; i++) {
		float r;
		r = erand48(dr48_buffer);
#if NOISE
		n[i] = 1e-4*r;
#else
		n[i] = 0.0*r;
#endif
	}
}

float Loss(float *x, int p, unsigned short *dr48_buffer)
{
	float n[p];
	float sum = 0.0;

	Noise(x, n, p, dr48_buffer);

	inc_nfc();
	sum = fitfun(x, p, NULL, NULL);
	sum += n[0];

#if VERBOSE
	printf("loss(%lf,%lf,%lf,%lf) = %lf\n", x[0], x[1], x[2], x[3], sum);
#endif
	return sum;
}


float *ghatinput;	// accumulated gradient vector
float **Hhatinput;	// accumulated hessian matrix
pthread_mutex_t hat_m = PTHREAD_MUTEX_INITIALIZER;	// lock for protecting multithreaded updates
//int gH_avg;		// random vectors generated (per task)

void do_reduce()
{
	if (torc_num_nodes() == 1) return;
#if VERBOSE
	int i;
	for (i = 0; i < data.Nth; i++) 
		printf("ghatinput[%d] = %f\n", i, ghatinput[i]);
#endif

	if (torc_node_id() == 0) {
		MPI_Reduce(MPI_IN_PLACE, ghatinput, data.Nth, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(MPI_IN_PLACE, Hhatinput[0], data.Nth*data.Nth, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#if VERBOSE
		printf("after\n");
		for (i = 0; i < data.Nth; i++) 
			printf("ghatinput[%d] = %f\n", i, ghatinput[i]);
#endif
	}
	else {
		MPI_Reduce(ghatinput, ghatinput, data.Nth, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(Hhatinput[0], Hhatinput[0], data.Nth*data.Nth, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
}

void do_work(int *pid)
{
	int p = data.Nth; 
	float c = data.uh; // same diffstep for every direction //1e-5;
//	int gH_avg = 20;	//%no. of random vectors generated 

	float *theta = data.theta;
	int gH_avg = data.reps / torc_num_workers();

//%%% GENERATION OF AVERAGED GRADIENT AND HESSIAN %%% %%% (NO AVERAGING IF gH_avg=1) %%%
	int id = *pid;

	unsigned short dr48_buffer[3];

	if (data.seed == 0) {
		dr48_buffer[0] = 0;
		dr48_buffer[1] = 0;
		dr48_buffer[2] = time(0);
	}
//	srand48_r(id, &dr48_buffer);
//	srand48_r((id+1)*time(0), &dr48_buffer);
//	srand48_r((id+1)*data.seed, &dr48_buffer);
//	dr48_buffer[2] = (id);
//	dr48_buffer[2] = (id+1)*time(0);
	dr48_buffer[2] = (id+1)*data.seed;

	int i, j;
	int m;
	for (m=1; m <=gH_avg; m++) {
		float delta[p];

		for (i = 0; i < p; i++) {
			float r;
			r = erand48(dr48_buffer);
			delta[i] = 2*round(r)-1;
		}

		float thetaplus[p], thetaminus[p];
		for (i = 0; i < p; i++) {
			thetaplus[i] = theta[i] + c*delta[i];		
			thetaminus[i] = theta[i] -c*delta[i];
		}

		float yplus  = Loss(thetaplus, p, dr48_buffer);
		float yminus = Loss(thetaminus, p, dr48_buffer);
            
		float ghat[p];
		for (i = 0; i < p; i++) 
			ghat[i]=(yplus-yminus)/(2*c*delta[i]);
            

		// update rank-level information
		pthread_mutex_lock(&hat_m);
		for (i = 0; i < p; i++) {
			ghatinput[i] += ghat[i];
		}
		pthread_mutex_unlock(&hat_m);

		//%%% GENERATE THE HESSIAN UPDATE %%% 
		float deltatilda[p];

		for (i = 0; i < p; i++) {
			float r;
			r = erand48(dr48_buffer);
			deltatilda[i]=2*round(r)-1;
		}

		float thetaplustilda[p], thetaminustilda[p];
		float thetaplusminustilda[p], thetaminusplustilda[p];

		for (i = 0; i < p; i++) {
			thetaplustilda[i] = thetaplus[i] + c*deltatilda[i];		
			thetaminustilda[i] = thetaminus[i] -c*deltatilda[i];

			thetaplusminustilda[i] = thetaplus[i] - c*deltatilda[i];		
			thetaminusplustilda[i] = thetaminus[i] + c*deltatilda[i];
		}

		float yplustilda=Loss(thetaplustilda, p, dr48_buffer);
		float yminustilda=Loss(thetaminustilda, p, dr48_buffer);
		float yplusminustilda=Loss(thetaplusminustilda, p, dr48_buffer);
		float yminusplustilda=Loss(thetaminusplustilda, p, dr48_buffer);
            
		float ghatplus[p], ghatminus[p];
		for (i = 0; i < p; i++) {
			ghatplus[i]  = (yplustilda-yplusminustilda)/(2*c*deltatilda[i]); 
			ghatminus[i] = (yminusplustilda-yminustilda)/(2*c*deltatilda[i]);
		}
            
		//%%% STATEMENT PROVIDING AN AVERAGE OF SP GRAD. APPROXS.PER ITERATION %%%
           
		float deltaghat[p];
		for (i = 0; i < p; i++) {
			deltaghat[i]=ghatplus[i]-ghatminus[i];
		}

		float Hhat[p][p];
		float HhatT[p][p];

		//for i=1:p
		//	 Hhat(:,i)=deltaghat(i)./(c*delta);
		//end
		for (j=0; j<p; j++) {
			for (i=0; i<p; i++) {
				Hhat[i][j]=deltaghat[j]/(2*c*delta[i]);
				HhatT[j][i] = Hhat[i][j];
			}
		}

		//Hhat=.5*(Hhat+Hhat');
		for (i=0; i<p; i++) {
			for (j=0; j<p; j++) {
				Hhat[i][j]=0.5*(Hhat[i][j] + HhatT[i][j]);
			}
		}

		// update rank-level information
		pthread_mutex_lock(&hat_m);
		for (i=0; i<p; i++) {
			for (j=0; j<p; j++) {
				Hhatinput[i][j]+=Hhat[i][j];
			}
		}
		pthread_mutex_unlock(&hat_m);


	}
}

int main(int argc, char *argv[])
{
	data_init();

	data.Nth = data.Nth;

	//theta = allocate1D(p);	//theta=zeros(p,1);

	ghatinput = allocate1D(data.Nth);
	Hhatinput = allocate2D(data.Nth);

	int i,j;

	torc_register_task(do_work);
	torc_register_task(do_reduce);
	aux_init();

	torc_init(argc, argv, 0);

	fprint_matrix_1d(stdout, "theta", data.theta, data.Nth);

//%%% Loop for number of hessians %%%

	reset_nfc();
	int ntasks = torc_num_workers();
	for (i = 0; i < ntasks; i++) {
		torc_create(-1, do_work, 1,
				1, MPI_INT, CALL_BY_COP, &i);
	}
	torc_waitall();


	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, do_reduce, 0);
	}
	torc_waitall();

	//printf("after reduction\n");

	/* compute the average values */
	int gH_avg = data.reps / torc_num_workers();

	for (i=0; i<data.Nth; i++) {
		for (j=0; j<data.Nth; j++) {
			Hhatinput[i][j] /= (ntasks*gH_avg);
		}
		ghatinput[i] /= (ntasks*gH_avg);
	}

	/* store them */
	float *gbar = allocate1D(data.Nth);
	float **Hbar = allocate2D(data.Nth);

	for (i=0; i<data.Nth; i++) {
		for (j=0; j<data.Nth; j++) {
			Hbar[i][j] = Hhatinput[i][j];
		}
		gbar[i] = ghatinput[i];
	}


	//%%% THETA UPDATE BELOW USES GAUSSIAN ELIMINATION %%% %%% TO AVOID DIRECT COMPUTATION OF HESSIAN INVERSE %%%

	float **Hbarbar = allocate2D(data.Nth); //[PROBDIM][PROBDIM];
	float **sHbarbar = allocate2D(data.Nth);

//	Hbarbar=sqrtm(Hbar*Hbar+0.000001*eye(data.Nth)/1.0); 
	int l;
	for (i=0; i<data.Nth; i++) {
		for (j=0; j<data.Nth; j++) {
			float s=0;
			for (l=0; l<data.Nth; l++) {
				s+= Hbar[i][l]*Hbar[l][j];
			}
			Hbarbar[i][j] = s;
			if (i==j) Hbarbar[i][j] += 0.000001;
		}
	}
	if (sqrtm(Hbarbar, sHbarbar, data.Nth)) {
		return 1;
	}

//	fprint_matrix_2d(stdout, "HESSIAN", sHbarbar, data.Nth, data.Nth);

	FILE *fp = fopen("sHbarbar.txt", "w");
	fprint_matrix_2d(fp, "sHbarbar", sHbarbar, data.Nth, data.Nth);
	fclose(fp);


      	printf("GRADIENT VECTOR :\n");
        for (i = 0; i < data.Nth; i++) {
                printf("%15.8lf ", gbar[i]);
        }
	printf("\n"); fflush(0);
        printf("HESSIAN MATRIX :\n");
        for (i = 0; i < data.Nth; i++) {
                for (j = 0; j < data.Nth; j++)
                        printf("%15.8lf ", sHbarbar[i][j]);
                printf("\n");
        }
        printf("\n"); fflush(0);

	gsl_matrix *hessian_mat = gsl_matrix_alloc(data.Nth, data.Nth);
	for(i=0; i<data.Nth; i++){
		for(j=0; j<data.Nth; j++){
			gsl_matrix_set(hessian_mat, i, j, sHbarbar[i][j]);
		}
	}
	eigs(hessian_mat, data.Nth);

//	print(eig(Hbarbar));

	free(data.theta);
	free(Hbar);
	free(gbar);
	free(ghatinput);
	free(Hhatinput);


	get_nfc();
	printf("total function calls = %d\n", get_tfc());

	torc_finalize();
	return 0;
}


//http://yarchive.net/comp/sqrtm.html
int sqrtm(float **A, float **sA, int n)
{
	int i, j;


#if VERBOSE
	fprint_matrix_2d(stdout, "A", A, n, n);
#endif
	gsl_matrix *mA = gsl_matrix_alloc (n, n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			gsl_matrix_set(mA, i, j, A[i][j]);
		}
	}

#if 0
	printf("m =\n ");
	gsl_matrix_fprintf(stdout, mA, "%.4lf");
#endif


	gsl_vector *Eig = gsl_vector_alloc (n);
	gsl_matrix *V = gsl_matrix_alloc (n, n);

	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
	gsl_eigen_symmv (mA, Eig, V, w);
	gsl_eigen_symmv_free (w);
//	gsl_eigen_symmv_sort (Ev, V, GSL_EIGEN_SORT_ABS_ASC);

	gsl_matrix *D = gsl_matrix_alloc (n, n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			gsl_matrix_set(D, i, j, 0.0);
		}
	}

	for (i = 0; i < n; i++) {
		float eval_i = gsl_vector_get (Eig, i);
		gsl_matrix_set(D, i, i, sqrt(eval_i));
	}
  
#if VERBOSE
	printf("m =\n ");
	gsl_matrix_fprintf(stdout, mA, "%.4lf");
	printf("V =\n ");
	gsl_matrix_fprintf(stdout, V, "%.4lf");
	printf("D =\n ");
	gsl_matrix_fprintf(stdout, D, "%.4lf");
#endif

	gsl_matrix *S0 = gsl_matrix_alloc (n, n);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
			1.0, V, D,
			0.0, S0);

	gsl_matrix *S = gsl_matrix_alloc (n, n);

	gsl_blas_dgemm (CblasNoTrans, CblasTrans,
			1.0, S0, V,
			0.0, S);

#if VERBOSE
	printf("S =\n ");
	gsl_matrix_fprintf(stdout, S, "%.4lf");
#endif

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			float v = gsl_matrix_get(S, i, j);
			sA[i][j] = v;
		}
	}

	gsl_vector_free (Eig);
	gsl_matrix_free (V);
	gsl_matrix_free (D);
	gsl_matrix_free(S0);
	gsl_matrix_free(S);

	return 0;	// everything ok
}

float** allocate2D(int n)
{
	int i;
#if 0
	float **a = (float **)calloc(1, n * sizeof(float *) + (n * n * sizeof(float)));
	float *mem = (float *)(a + n);
#else
	float *mem = (float *)calloc(1, n*n*sizeof(float));
	float **a =  (float **) malloc(n*sizeof(float *));
#endif
	for (i = 0; i < n; i++) {
		a[i] = mem + (i*n); 
	}

	return a;
}

float* allocate1D(int n)
{
	float *a = (float *)calloc(1, n * sizeof(float));
	return a;
}

