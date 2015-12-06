#ifndef _MARKET_FIXED_POINT_
#define	_MARKET_FIXED_POINT_

typedef struct{
	void (*function)(double *in, void * p, double *out);
	void *params;
} ff_function;

void fixed_point(int n, ff_function *f, double *x0,
                 double xtol, int maxiter, double *out);

#endif
