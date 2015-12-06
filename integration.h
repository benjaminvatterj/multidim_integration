#ifndef _MARKET_INTEGRATION_
#define _MARKET_INTEGRATION_
#include <gsl/gsl_integration.h>

typedef struct{
	double (*function)(double *x, size_t n, void * p);
	void *params;
} nquad_function;

typedef void (*multi_integration_type)(int n, nquad_function *f, double *lb, double *ub, double epsabs,
           double epsrel, size_t limit, double *result, double *abserr);

void nquad(int n, nquad_function *f,  double *lb, double *ub,
           double epsabs, double epsrel, size_t limit,
           double * result, double *abserr);

void mcint(int n, nquad_function *f, double *lb, double *ub, double epsabs,
           double epsrel, size_t limit, double *result, double *abserr);


#endif
