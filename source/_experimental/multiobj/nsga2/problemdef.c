/* Test problem definitions */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <unistd.h>
# include "global.h"
# include "rand.h"

void fitfun(float /*const*/ *x, int N, void *output, int *winfo, float *result, int ncon, float *constraints);

void test_problem (float *xreal, float *xbin, int **gene, float *obj, float *constr)
{
	obj[0] = pow((xreal[0]+2.0),2.0)-10;
	obj[1] = pow((xreal[0]-2.0),2.0)+20;
	return;
}

void test_problem_v2 (float *x, int *pnx, float *obj, int *pno, float *constr, int *pncon, float *constr_violation, int info[4])
{
	int nx = *pnx;
	int no = *pno;
	int ncon = *pncon;

	//int popid = info[2];
	//int id = info[3];

	//usleep(10*1000);

	// Option #0 (original): everything here
	//obj[0] = pow((x[0]+2.0),2.0)-10;
	//obj[1] = pow((x[0]-2.0),2.0)+20;
	//constr[0] = obj[1] - 25;

	// Option #1: fitfun code moved to fitfun.c and constraints stay here
	fitfun(x, nx, NULL, info, obj, ncon, constr);
	//constr[0] = obj[1] - 25;

	// Option #2: fitfun also handles the constraints
	//fitfun(x, nx, NULL, info, obj, ncon, constr);

	// the following stays as it was
	if (ncon==0)
	{
		*constr_violation = 0.0;
	}
	else
	{
		*constr_violation = 0.0;

		int j;
		for (j=0; j<ncon; j++)
		{
			if (constr[j]<0.0)
			{
				*constr_violation += constr[j];
			}
		}
	}

        return;
}


void test_problem_v2_old (float *x, int *pnx, float *obj, int *pno, int info[4])
{
	int nx = *pnx;
	//int no = *pno;

	//int popid = info[2];
	//int id = info[3];

	//usleep(10*1000);
	//obj[0] = pow((x[0]+2.0),2.0)-10;
	//obj[1] = pow((x[0]-2.0),2.0)+20;
	fitfun(x, nx, NULL, info, obj, 0, NULL);

//	obj[0] = 4.0*pow(x[0],2.0)+4.0*pow(x[1],2.0);
//	obj[1] = pow((x[0]-5.0),2.0)+ pow((x[1]-5.0),2.0);

	return;
}
