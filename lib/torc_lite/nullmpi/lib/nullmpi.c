#include <stdio.h>
#include <sys/time.h>

void empty()
{
}

float MPI_Wtime()
{
	struct timeval t;      
	gettimeofday(&t, NULL);
	return (float)t.tv_sec + (float)t.tv_usec*1.0E-6;
}

