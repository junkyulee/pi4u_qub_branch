#ifndef _PNDLC_H_
#define _PNDLC_H_

void c_pndlhfa(float (*F)(float *, int *), const float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
                                const int *iprint, float *hes, const int *ld, int *noc, int *ierr);
void c_pndlga(float (*F)(float *, int *), const float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
                                const int *iprint, float *grad, int *noc, int *ierr);
void c_pndlghfa(float (*F)(float *, int *), const float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
                                const int *iprint, float *grad, float *hes, const int *ld, int *noc, int *ierr);
void c_pndl_init(void);
void c_pndl_finalize(void);
void c_pndlhga(void (*GRD)(float *, int *, float *), float *x, const int *n, const float *xl, const float *xu, const float *uh, const float *feps, const int *iord,
                                const int *iprint, float *hes, int *ld, int *noc, int *ierr);
void c_pndlja(void (*RSD)(float *, int *, int *, float *), float *x, const int *n, const int *m, const float *xl, const float *xu, const float *uh, const float *feps,
                                const int *iord, const int *iprint, float *fj, int *ld, int *noc, int *ierr);
#endif
