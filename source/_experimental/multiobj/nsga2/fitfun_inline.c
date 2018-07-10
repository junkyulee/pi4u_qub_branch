#include <math.h>
#include <unistd.h>

void fitfun(float *x, int nx, void *output, int *winfo, float *result, int ncon, float *constraints) // ny must be also available
{
	usleep(10*1000);

//	result[0] = pow((x[0]+2.0),2.0)-10;
//	result[1] = pow((x[0]-2.0),2.0)+20;

	float y = 0;

	int i;
	for (i = 0; i < nx; i++) y += x[i];

	result[0] = pow((y+2.0),2.0)-10;
	result[1] = pow((y-2.0),2.0)+20;
}

