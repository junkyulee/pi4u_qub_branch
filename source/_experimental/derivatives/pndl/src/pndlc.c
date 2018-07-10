#ifndef _PNDLC_H_
#define _PNDLC_H_

#include "pndl_config.h"

void F77_FUNC(pndlhfa,PNDLHFA)(float (*F)(float *, int *), const float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
                                const int *iprint, float *hes, const int *ld, int *noc, int *ierr);
void F77_FUNC(pndlga,PNDLGA)(float (*F)(float *, int *), const float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
                                const int *iprint, float *grad, int *noc, int *ierr);
void F77_FUNC(pndlghfa,PNDLGHFA)(float (*F)(float *, int *), const float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
                                const int *iprint, float *grad, float *hes, const int *ld, int *noc, int *ierr);
void F77_FUNC(pndl_init,PNDL_INIT)(void);
void F77_FUNC(pndl_finalize,PNDL_FINALIZE)(void);
void F77_FUNC(pndlhga,PNDLHGA)(void (*GRD)(float *, int *, float *), float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
				const int *iprint, float *hes, int *ld, int *noc, int *ierr);
void F77_FUNC(pndlja,PNDLJA)(void (*RSD)(float *, int *, int *, float *), float *x, const int *n, const int *m, const float *xl, const float *xu, const float *uh, const float *feps,
				const int *iord, const int *iprint, float *fj, int *ld, int *noc, int *ierr);


void c_pndlhfa(float (*F)(float *, int *), const float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
                                const int *iprint, float *hes, const int *ld, int *noc, int *ierr)
{
	F77_FUNC(pndlhfa,PNDLHFA)(F, x, n, xl, xu, uh, feps, iord, iprint, hes, ld, noc, ierr);
}

void c_pndlga(float (*F)(float *, int *), const float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
                                const int *iprint, float *grad, int *noc, int *ierr)
{
	F77_FUNC(pndlga,PNDLGA)(F, x, n, xl, xu, uh, feps, iord, iprint, grad, noc, ierr);
}

void c_pndlghfa(float (*F)(float *, int *), const float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
                                const int *iprint, float *grad, float *hes, const int *ld, int *noc, int *ierr)
{
	F77_FUNC(pndlghfa,PNDLGHFA)(F, x, n, xl, xu, uh, feps, iord, iprint, grad, hes, ld, noc, ierr);
}

void c_pndl_init(void)
{
	F77_FUNC(pndl_init,PNDL_INIT)();
}

void c_pndl_finalize(void)
{
	F77_FUNC(pndl_finalize,PNDL_FINALIZE)();
}

void c_pndlhga(void (*GRD)(float *, int *, float *), float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
				const int *iprint, float *hes, int *ld, int *noc, int *ierr)
{
	F77_FUNC(pndlhga,PNDLHGA)(GRD, x, n, xl, xu, uh, feps, iord, iprint, hes, ld, noc, ierr);
}

void c_pndlja(void (*RSD)(float *, int *, int *, float *), float *x, const int *n, const int *m, const float *xl, const float *xu, const float *uh, const float *feps,
				const int *iord, const int *iprint, float *fj, int *ld, int *noc, int *ierr)
{
	F77_FUNC(pndlja,PNDLJA)(RSD, x, n, m, xl, xu, uh, feps, iord, iprint, fj, ld, noc, ierr);
}

#endif
